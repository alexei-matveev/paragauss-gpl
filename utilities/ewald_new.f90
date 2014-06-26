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
     module type_module
       integer, parameter :: r8_kind=selected_real_kind(15)
       integer, parameter :: r4_kind=selected_real_kind(6)
       integer, parameter :: i4_kind=selected_int_kind(9)
     end module type_module
!
     module error_module
       contains
       subroutine error_handler(message)
         implicit none
         character(len=*), intent(in) :: message
         write(*,'(10X,60(1h*)/15x,a50/30X,"Terminated"/10X,60(1h*))') &
              adjustl(trim(message))
         stop
       end subroutine error_handler
     end module error_module
!
! 
!
     module ewald_module
       use type_module
       use error_module
       implicit none

       type screep_point_type
          real(kind=r8_kind), dimension(3) :: rc
          real(kind=r8_kind)               :: V_dif,q,square
       end type screep_point_type

       type atom_type
	  real(kind=r8_kind), dimension(3) :: s    ! shell  coordinate
	  real(kind=r8_kind), dimension(3) :: c    ! core  coordinate
          real(kind=r8_kind), dimension(3) :: rc   ! cartesian coordinate
          real(kind=r8_kind), dimension(3) :: rl   ! lattice   coordinate
          real(kind=r8_kind)               :: q    ! charge
          integer(kind=i4_kind)            :: qt   ! charge type
	  real(kind=r8_kind)               :: an   ! real for charge type
          character(len=6)                 :: name ! atom symbol
       end type atom_type

       type cluster_type
          real(kind=r8_kind), dimension(3) :: rc   ! cartesian coordinate
          real(kind=r8_kind), dimension(3) :: rl   ! lattice   coordinate
          real(kind=r8_kind)               :: V_pc ! PC-array potential on ato          
          real(kind=r8_kind)               :: V_ew ! madelung potential on atom
          character(len=6)                 :: name ! atom symbol
       end type cluster_type

       type lattice_type
          real(kind=r8_kind),    pointer, dimension(:,:) :: v_cartesian
          real(kind=r8_kind),    pointer, dimension(:)   :: v_mod
          integer(kind=i4_kind), pointer, dimension(:,:) :: v_vector
          integer(kind=i4_kind)                          :: num_vector
          integer(kind=i4_kind), pointer, dimension(:)   :: num_in_sphere
          integer(kind=i4_kind)                          :: n_spheres
       end type lattice_type

       type point_charge
          real(kind=r8_kind), dimension(3) :: rc
          real(kind=r8_kind)               :: q
          character(len=6)                 :: name
       end type point_charge
       
       type(screep_point_type), pointer, dimension(:)     :: screep_point
       type(atom_type),         pointer, dimension(:)     :: atom
       type(cluster_type),      pointer, dimension(:)     :: cluster
       type(point_charge),      pointer, dimension(:)     :: pc
       real(kind=r8_kind),  allocatable, dimension(:,:)   :: vertex,rl
       real(kind=r8_kind),  allocatable, dimension(:)     :: vlm
       real(kind=r8_kind)                  :: R_screep, V_cell, V_BZ, eta
       real(kind=r8_kind), dimension(3)    :: Rc_screep, &
            Rc_screepAU,Rc_screepANGS,direction
       real(kind=r8_kind), dimension(3,3)  :: a, b
       real(kind=r8_kind), parameter       :: pi=3.14159265358979_r8_kind
       integer(kind=i4_kind)               :: n_triangle, n_atoms_in_cell, n_cluster_atoms
       integer(kind=i4_kind), parameter    :: input_unit=5, output_unit=6 ,gxc_unit=7,gx_unit=8
       integer(kind=i4_kind)               :: N_ewald_pc,n_screep_points, n_div
       integer(kind=i4_kind)               :: i_scheme,lmax
       character(len=100)                  :: work_dir, input_file, output_file
       character(len=256)                  :: libdir
       character(len=100)                  :: screep_procedure
       logical                             :: axes_c3, axes_c4, input_file_exist
       logical                             :: symmetry_c3v=.false.

       type(lattice_type)                  :: lattice
       real(kind=r8_kind)                  :: R_max_pc, R_max_fragment
       logical, parameter                  :: print_vectors=.false. !print_vectors=.true.
       logical, parameter                  :: not_print_vectors=.false.
       logical, parameter                  :: full_print=.false. !full_print=.true.
       logical                             :: make_crystal_fragment, make_ewald_array ! tasks
       logical:: gx_file_exist
       logical:: gxcv_file_exist

	logical, dimension(:,:), allocatable:: type_number

     contains

       subroutine ewald

         implicit none         
         character(len=256)                  :: gxcv_file
         character(len=256)                   :: gx_cluster
         real(kind=r8_kind), pointer,dimension(:)        :: V_ew    ,qp
         real(kind=r8_kind), allocatable, dimension(:,:) :: A,c
         real(kind=r8_kind)                              :: R,zero=0.0_r8_kind,one=1.0_r8_kind
         real(kind=r8_kind)                              :: cl_V
         integer(kind=i4_kind)                           :: i_atom,i_point,alloc_stat,i,j
         integer(kind=i4_kind)                           :: nlm,m,l,ns,nl
         character(len=100)                              :: format

!         call getenv_("PWD",work_dir) undescored getenv_(...) for Linux
         call getenv("PWD",work_dir)
         call getenv("TTFSLIBS",libdir)
         input_file= trim(work_dir)//"/"//"epe.input"
         output_file=trim(work_dir)//"/"//"ewald.output"
	 gxcv_file=trim(work_dir)//"/"//"cellvec"
         nl=len_trim(work_dir)
         ns=index(trim(work_dir),"/",back=.true.)
         gx_cluster=trim(work_dir)//"/"//"gx."//work_dir(ns+1:nl)

         inquire(file=input_file, exist=input_file_exist)
         inquire(file=gxcv_file, exist=  gxcv_file_exist)
	 if(gxcv_file_exist)  then
            open(gxc_unit,file=trim(gxcv_file),status="old")
         endif
         if(.not.input_file_exist) then
            write(*,*) "Making template file : ", input_file
            call error_handler&
                 ("Screep_input : file <ewald_input> not exist !!!")
         end if

         inquire(file=gx_cluster, exist=gx_file_exist)
         if(.not.gx_file_exist)  then
            gx_cluster=trim(work_dir)//"/"//"gx"
         end if
         inquire(file=gx_cluster, exist=gx_file_exist)
         if(.not.gx_file_exist)  then
            gx_cluster=trim(work_dir)//"/"//"gxfile"
         end if
         inquire(file=gx_cluster, exist=gx_file_exist)
         if(gx_file_exist)  then
            open(gx_unit,file=trim(gx_cluster),status="old")
         end if

         open(input_unit,  file=trim(input_file), status="unknown")
         open(output_unit, file=trim(output_file),status="unknown")
         write(output_unit,*) "Input directory   : ", trim(work_dir)
         write(output_unit,*) "Input  file       : ", trim(input_file)
         write(output_unit,*) "Output file       : ", trim(output_file)
         write(output_unit,*) "FF directory      : ", trim(libdir)
         if(gxcv_file_exist) &
              write(output_unit,*) "External unit cell: ",trim(gxcv_file)
         if(gx_file_exist) &
              write(output_unit,*) "QM cluster        : ",trim(gx_cluster)

         call ewald_input_new
         call pc_generate 
         if(make_crystal_fragment) stop ! terminate work after the first stage is performed
         allocate(rl(n_cluster_atoms,3), stat=alloc_stat)
         if(alloc_stat.ne.0) &
              call error_handler("Ewald error [1]: Allocation failed")
         do i_atom=1,n_cluster_atoms
            rl(i_atom,:)=cluster(i_atom)%rl(:)
         end do
         V_ew => cluster(:)%V_ew
         call Madelung(rl,V_ew,n_cluster_atoms)

         deallocate(rl, stat=alloc_stat)
         if(alloc_stat.ne.0) &
              call error_handler("Ewald error [2]: Deallocation failed")
         do i_atom=1,n_cluster_atoms
            write(output_unit,*) i_atom, "V=  ",cluster(i_atom)%V_ew 
         end do
         
         select case(screep_procedure)
            case("screep_pc")
               call screep_generate
            case("screep_multipole")
               call genmesh(Rc_screep,r_screep,i_scheme)
            case("screep_multipole_pc")
         end select

         allocate(rl(n_screep_points,3), stat=alloc_stat)
         if(alloc_stat.ne.0) &
              call error_handler("Ewald error [3]: Allocation failed")
         do i_point=1,n_screep_points
            rl(i_point,:)=matmul(screep_point(i_point)%rc,b)/(2*pi)
         end do
         V_ew => screep_point(:)%V_dif
         call Madelung(rl,V_ew,n_screep_points)
         deallocate(rl, stat=alloc_stat)
         if(alloc_stat.ne.0) &
              call error_handler("Ewald error [4]: Deallocation failed")

         do i_point=1,n_screep_points
            do i_atom=1,N_ewald_pc
               R=sqrt( dot_product(screep_point(i_point)%rc-pc(i_atom)%rc,&
                                   screep_point(i_point)%rc-pc(i_atom)%rc ) )
               if(R<=1.d-10) R=zero
               if(R /= zero) then
                  screep_point(i_point)%V_dif=screep_point(i_point)%V_dif &
                                                            - pc(i_atom)%q/R
               end if
            end do
         end do

         select case(screep_procedure)
         case("screep_pc")
            allocate(A(n_screep_points, n_screep_points), stat=alloc_stat)
            if(alloc_stat.ne.0) &
                 call error_handler("Ewald error [5]: Deallocation failed")
            do i=1, n_screep_points
               A(i,i)=1.07d0 * sqrt(4*acos(-one)/screep_point(i)%square)
               do j=i+1, n_screep_points
                  R=sqrt(dot_product(screep_point(i)%rc-screep_point(j)%rc, &
                                     screep_point(i)%rc-screep_point(j)%rc) )
                  A(i,j)=one/R
                  A(j,i)=A(i,j)
               end do
            end do
            call inverse_matrix(A,n_screep_points)               
            V_ew => screep_point(:)%V_dif
            screep_point(:)%q = matmul(A,V_ew)

            write(output_unit,"(/10X, 'Total charge on sphere :', f9.5)") &
                    sum(screep_point(:)%q)
            do i=1,n_screep_points
               write(output_unit,*) i, 'q=', screep_point(i)%q
            end do
            write(output_unit,*) ' '
            write(output_unit,*) '                 V_ew         V_cl'
            do i_atom=1,n_cluster_atoms
               cl_V=cluster(i_atom)%V_pc
               do j=1,n_screep_points
                  R=sqrt(dot_product(screep_point(j)%rc-cluster(i_atom)%rc, &
                                     screep_point(j)%rc-cluster(i_atom)%rc  ) )
                  if(R <= 1.d-10) R=zero
                  if(R /= zero) cl_V=cl_V + screep_point(j)%q/R
               end do
               format="(10x,i3,3x,3(f10.7,3x))"
               write(output_unit,format) i_atom,cluster(i_atom)%V_ew, cl_V, &
                                                abs(cluster(i_atom)%V_ew-cl_V)
            end do
!++++++++++++++++++++++++++++ Store embedded cluster ++++++++++++++++++++++++++
            open(20, file="ew_pc.xyz", status="unknown")
            open(30, file="ewald.pcr", status="unknown")
            write(20,*) n_screep_points+n_ewald_pc
	    write(30,*) n_screep_points+n_ewald_pc
            write(20,*) "comment line"
            i_at: do i_atom=1,n_ewald_pc
               do i=1,n_cluster_atoms
                  if( all(abs(pc(i_atom)%rc-cluster(i)%rc)<=1.d-10) ) then
                     write(20,"(' Si',4x,3(f12.7,3x),15x,3x,i4)") &
                             0.529177_r8_kind* (cluster(i)%rc(:)+Rc_screepAU), i
		    write(30,"(4f15.8,i3,4i2,6x,f4.2)") &
                    pc(i_atom)%rc(:)+Rc_screepAU, pc(i_atom)%q,1,0,0,1,0, 0.0 
                     cycle i_at
                  end if
               end do
               write(20,"(a4,4x,3(f12.7,3x),5x,f10.7,3x,i4)") &
                    pc(i_atom)%name, &
                    0.529177_r8_kind*( pc(i_atom)%rc(:)+Rc_screepAU), pc(i_atom)%q, i_atom   
               write(30,"(4f15.8,i3,4i2,6x,f4.2)") &
		    pc(i_atom)%rc(:)+Rc_screepAU, pc(i_atom)%q,1,0,0,1,0, 0.0 
            end do i_at

            do i=1,n_screep_points
               write(20,"('H',3(5x,f12.7))" ) &
                    (RC_screepAU+screep_point(i)%rc(:))*0.529177_r8_kind
	      write(30,"(4f15.8,i3,4i2,6x,f4.2)") &
	      screep_point(i)%rc(:)+Rc_screepAU, &
              screep_point(i)%q,1,0,0,1,0, 0.0
            end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case("screep_multipole")
            allocate(c(lmax+1,2*lmax+1), vlm((lmax+1)**2),stat=alloc_stat)
            if(alloc_stat.ne.0) &
                 call error_handler("Ewald error [6]: Allocation failed")

            V_ew => screep_point(:)%V_dif
            call project(V_ew,lmax,Rc_screep,c)
            do i_atom=1,n_cluster_atoms
               call rest(cluster(i_atom)%rc(1),cluster(i_atom)%rc(2), &
                         cluster(i_atom)%rc(3),r_screep,Rc_screep,lmax,c,cl_V)
               cl_V=cl_V + cluster(i_atom)%V_pc
               format="(10x,i3,3x,3(f10.7,3x))"
               write(output_unit,format) i_atom, cluster(i_atom)%V_ew, cl_V, &
                                                 abs(cluster(i_atom)%V_ew-cl_V)
            end do
            deallocate(c,vlm, stat=alloc_stat)
            if(alloc_stat.ne.0) &
                 call error_handler("Ewald error [7]: Deallocation failed")
         end select
         nullify(V_ew)

       end subroutine ewald
       
!#####################################################
       subroutine ewald_input_new

         implicit none

         integer(i4_kind) :: line_number
         character*80 :: current_block
         logical :: read_block,end_block,comment_line,empty_line
         character*10 :: files(2)
         integer(i4_kind), parameter :: epeinp=1, fflib=2
         character*30 :: ff_name
         integer(i4_kind) :: ff_unit=9
         logical :: l_ff

         real(kind=r8_kind), parameter ::zero=0.0_r8_kind,two=2.0_r8_kind, three=3.0_r8_kind, &
              one=1.0_r8_kind
         real(kind=r8_kind)::an
         real(kind=r8_kind), allocatable, dimension(:) :: tmp_r
         integer(kind=i4_kind)                         :: alloc_stat,ioerr
         character(len=100)                            :: format
         character(len=100)                            :: axes,units, cluster_key,pc_key, task_key
         real(kind=r8_kind) :: u(3,3), axis(3),pos(3),pos2(3),ang(3)
         real(kind=r8_kind), allocatable :: r_buf(:,:)
         integer(kind=i4_kind), allocatable :: perm(:,:)
         logical :: new_reg,read_at
         integer(kind=i4_kind) :: max_type_ions
         character(len=6) :: name_of_type(100)

         files(1)="epeinput"; files(2)="forcefield"

         make_crystal_fragment=.false.
         make_ewald_array=.true.
         !Default values
         axes_c3 =.false.
         axes_c4 =.false.
         Rc_screep=(/0.d0,0.d0,0.d0/)
         r_screep=5.0_r8_kind
         screep_procedure="screep_pc"
         direction=(/0.d0,0.d0,1.d0/)
         n_screep_points=128
         R_max_pc=15.0_r8_kind
         eta=0.1_r8_kind

         if(gxcv_file_exist) call read_atom_names

         call read_cell_vectors

         call read_coordinates

         call read_screep_sphere

         call read_pc_array

         call read_cluster

         call read_ewald_parameter

       contains
         !==============================================================
         subroutine read_atom_names()

           character*80 :: buffer
           character*10 :: buffer1
           integer(i4_kind) :: i,j

           new_reg=.true.

           read_block=find_block(input_unit,"&ATOM_NAMES",epeinp)
           if(.not.read_block) then
              read_at=.false.
           else
              read_at=.true.
              i=0; j=0
              do
                 call read_block_line(input_unit,buffer,end_block,comment_line,empty_line,epeinp)
                 if(empty_line) cycle
                 if(comment_line) cycle
                 if(end_block) exit
                 i=i+1
                 if(i==1) then
                    !reading in number of atoms in the unit cell
                    call upcase(buffer)
                    if(check_string(buffer,"N_ATOMS")) then
                       read(buffer,*,err=100)  buffer1, n_atoms_in_cell
                       if(check_string(buffer,"OLD")) then
                          new_reg=.false.
                       else if(check_string(buffer,"NEW")) then
                          new_reg=.true.
                       end if
                    else
                       call error_handler("Read_input: The first line in &ATOM_NAMES"// &
                            " block has to be number of atoms (N_ATOMS  NNN)")
                    end if

                    format="(/10X, 'Number of atoms in cell : ', i4)"
                    write(output_unit,format) n_atoms_in_cell
                    allocate(atom(n_atoms_in_cell), tmp_r(3),stat=alloc_stat)
                    if(alloc_stat.ne.0) &
                         call error_handler("Read_input: Allocation failed [33]")
                    atom(:)%s(1)=0.0_r8_kind
                    atom(:)%s(2)=0.0_r8_kind
                    atom(:)%s(3)=0.0_r8_kind
                    allocate(type_number(n_atoms_in_cell,n_atoms_in_cell),stat=alloc_stat)
                    type_number=.false.
                    if(alloc_stat.ne.0) &
                         call error_handler("Read_input: Allocation failed [33a]")
                 else if(i==2) then
                    call upcase(buffer)
                    if(check_string(buffer,"N_ATOM_TYPES")) then
                       read(buffer,*,err=100)  buffer1, max_type_ions
                    else
                       call error_handler("Read_input: The second line in &ATOM_NAMES"// &
                            " block has to be number of atomic types (N_ATOM_TYPES  NNN)")
                    end if
                 else
                    j=j+1
                    read(buffer,*,err=100)  name_of_type(j)
                 end if
              end do
           end if

           return
           
100        call wrong_line(epeinp)

         end subroutine read_atom_names
         !==============================================================
         subroutine read_cell_vectors

           character*80 :: buffer
           integer(i4_kind) :: i,j,k,i_vector

           if(read_at) goto 10

           read_block=find_block(input_unit,"&CELL_VECTORS",epeinp)
           if(.not.read_block) then
              call error_handler("Read_input: No cell vectors found")
           else
              j=0
              do 
                 call read_block_line(input_unit,buffer,end_block,comment_line,empty_line,epeinp)
                 if(empty_line) cycle
                 if(comment_line) cycle
                 if(end_block) exit
                 j=j+1
                 if(j > 3) exit
                 read(buffer,*,err=100) (a(i,j),i=1,3)
              end do
              if(j < 3) call error_handler("Read_input: Number of cell vectors < 3")
           end if

10         if(gxcv_file_exist) then
              do j=1,3
                 read(gxc_unit,*,err=200)(a(i,j),i=1,3)
              enddo
           end if

           format="(10X,'a_',i1,3(f15.10,3x))"
           write (output_unit,*) '=============================================='
           write(output_unit,*) 'Cell Vectors'
           if(gxcv_file_exist) write(output_unit,*) 'Vectors are taken from cellvec file'
           do i_vector=1,3
              write(output_unit,format) i_vector, a(:,i_vector)
           end do
           write(output_unit,*) '----------------------------------------------'

           a=a/0.529177_r8_kind

           V_cell=dot_product(a(:,1),cross_product(a(:,2),a(:,3)))
           if(V_cell < 0.0) then
              format="(/10X, 'WARNING : (a_1,[a_2,a_3]) Negative (May be redefine a_i order)' )"
              write(output_unit,format)
           end if
           !Reciprocal lattice vectors
           do i=1,3
              j=i+1
              if(j.gt.3) j=j-3
              k=j+1
              if(k.gt.3) k=k-3
              b(:,i)=2*pi*cross_product(a(:,j),a(:,k) ) / V_cell
           end do
           format="(/10X, 'Unit cell volume : ',f20.12,' (a.u)**3' )"
           write(output_unit,format) V_cell
           V_cell=abs( V_cell )

           format="(/10X, 'Reciprocal Lattice vectors (x y z)')"
           write(output_unit,format)
           format="(10X,'b_',i1,3(f15.10,3x))"
           do i_vector=1,3
              write(output_unit,format) i_vector, b(:,i_vector)
           end do
           V_BZ=8*pi**3/V_cell
           format="(/10X, '    V_BZ  volume : ',f20.12,' (a.u)**3' )"
           write(output_unit,format) V_BZ

           return

100        call wrong_line(epeinp)
200        call error_handler("Read_input: Wrong CELLVEC file")

         end subroutine read_cell_vectors
         !==============================================================
         subroutine read_coordinates

           character*80 :: buffer
           character*10 :: buffer1
           real(r8_kind) :: r_sh(3),ann
           integer(i4_kind) :: i,j,i_atom
           character*6 :: atname
           real(r8_kind) :: q_core_epe,q_shell_epe

           if(read_at) goto 20

           read_block=find_block(input_unit,"&COORDINATES",epeinp)
           if(.not.read_block) then
              call error_handler("Read_input: No cell coordinates found")
           else
              i=0; j=0
              do
                 call read_block_line(input_unit,buffer,end_block,comment_line,empty_line,epeinp)
                 if(empty_line) cycle
                 if(comment_line) cycle
                 if(end_block) exit
                 i=i+1
                 if(i==1) then
                    !reading in number of atoms in the unit cell
                    call upcase(buffer)
                    if(check_string(buffer,"N_ATOMS")) then
                       read(buffer,*,err=100)  buffer1, n_atoms_in_cell
                       if(check_string(buffer,"OLD")) then
                          new_reg=.false.
                       else if(check_string(buffer,"NEW")) then
                          new_reg=.true.
                       end if
                    else
                       call error_handler("Read_input: The first line in &COORINATES"// &
                            " block has to be number of atoms (N_ATOMS  NNN)")
                    end if

                    format="(/10X, 'Number of atoms in cell : ', i4)"
                    write(output_unit,format) n_atoms_in_cell
                    allocate(atom(n_atoms_in_cell), tmp_r(3),stat=alloc_stat)
                    if(alloc_stat.ne.0) &
                         call error_handler("Read_input: Allocation failed [3]")
                    atom(:)%s(1)=0.0_r8_kind
                    atom(:)%s(2)=0.0_r8_kind
                    atom(:)%s(3)=0.0_r8_kind
                    allocate(type_number(n_atoms_in_cell,n_atoms_in_cell),stat=alloc_stat)
                    type_number=.false.
                    if(alloc_stat.ne.0) &
                         call error_handler("Read_input: Allocation failed [3a]")
                 else
                    !reading atomic coordinates
                    j=j+1
                    read(buffer,*,err=100)  atom(j)%name,tmp_r
                    atom(j)%q=0.0_r8_kind
                    atom(j)%rc=tmp_r
                 end if
              end do
              if(j /= n_atoms_in_cell) call error_handler("Read_input: Check number of atoms!!!")
           end if

20         if(gxcv_file_exist) then
              write(output_unit,*) 'Read coordinates of cell from cellvec file'
              do i_atom=1,n_atoms_in_cell
                 read(gxc_unit,*,err=300) atom(i_atom)%an, tmp_r
                 atom(i_atom)%qt=int(atom(i_atom)%an)
                 type_number(atom(i_atom)%qt,i_atom)=.true.
                 atom(i_atom)%c=tmp_r

                 if(read_at) then
                    atom(i_atom)%q=0.0_r8_kind
                    atom(i_atom)%name=name_of_type(int(atom(i_atom)%an))
                 end if
              end do
              do
                 read(gxc_unit,*,end=22,err=300) an,tmp_r
                 if(an.eq. zero) exit
                 !                ian=an
                 do i_atom=1,n_atoms_in_cell
                    if(type_number(int(an),i_atom)) then
                       atom(i_atom)%s=tmp_r
                       type_number(int(an),i_atom) = .false.
                       exit
                    endif
                 enddo
              enddo
22            continue
              do i_atom=1,n_atoms_in_cell
                 if(dot_product(atom(i_atom)%s,atom(i_atom)%s).ne. zero) then
                    if(new_reg) then 
                       atom(i_atom)%rc=atom(i_atom)%c
                    else
                       atom(i_atom)%rc=(atom(i_atom)%c+atom(i_atom)%s)/2.0_r8_kind
                    end if
                 else
                    atom(i_atom)%rc=atom(i_atom)%c
                 endif
              enddo
           end if

           !Now define charges from FF file
           read_block=find_block(input_unit,"&FORCE_FIELD",epeinp)
           if(.not.read_block) then
              call error_handler("Read_input: No force field found")
           else
              i=0
              do
                 call read_block_line(input_unit,buffer,end_block,comment_line,empty_line,epeinp)
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
              open(ff_unit,file=trim(libdir)//"/"//trim(ff_name),  &
                   form='formatted', status='unknown')
           else
              call error_handler("Read_input: file "//trim(libdir)//"/"//trim(ff_name)//" not found")
           end if

           read_block=find_block(ff_unit,"&ATOM_PROPERTIES",fflib)
           if(.not.read_block) then
              call error_handler("Read_fflib: No atomic properties found")
           else
              i=0
              do
                 call read_block_line(ff_unit,buffer,end_block,comment_line,empty_line,fflib)
                 if(empty_line) cycle
                 if(comment_line) cycle
                 if(end_block) exit
                 i=i+1
                 read(buffer,*,err=100) atname, q_core_epe, q_shell_epe

                 do j=1,n_atoms_in_cell
                    if(trim(atname)==trim(atom(j)%name)) atom(j)%q=q_core_epe+q_shell_epe
                 end do
              end do
           end if

           write(output_unit,*) "============================================================"
           format="(i4,3x,i4, 3f20.7, f12.7)"
           do i_atom=1,n_atoms_in_cell
              write(output_unit,format) i_atom,atom(i_atom)%qt, atom(i_atom)%rc,atom(i_atom)%q 
              atom(i_atom)%rc=atom(i_atom)%rc/0.529177_r8_kind
              atom(i_atom)%rl= &
                   matmul(atom(i_atom)%rc,b)/(2.0_r8_kind*pi)
           end do
           write(output_unit,*) "============================================================"

           return

100        call wrong_line(epeinp)
200        call wrong_line(fflib)
300        call error_handler("Read_input: Wrong CELLVEC file")

         end subroutine read_coordinates
         !==============================================================
         subroutine read_screep_sphere

           character*80 :: buffer,buffer1
           integer(i4_kind) :: i_atom

           axes="C3"
           axes_c3=.true.
           n_triangle=8

           read_block=find_block(input_unit,"&SCREEP_DATA",epeinp)
           if(.not.read_block) then
!              call error_handler("Read_input: No screep data found")
              goto 10
           else
              do 
                 call read_block_line(input_unit,buffer,end_block,comment_line,empty_line,epeinp)
                 if(empty_line) cycle
                 if(comment_line) cycle
                 if(end_block) exit

                 call upcase(buffer)
                 if(check_string(buffer,"CENTER")) then
                    read(buffer,*,err=100) buffer1, Rc_screep
                 elseif(check_string(buffer,"RADIUS")) then
                    read(buffer,*,err=100) buffer1, r_screep
                 elseif(check_string(buffer,"C3V")) then
                    axes="C3V"
                    axes_c3=.true.
                    axes_c4=.false.
                    symmetry_c3v=.true.
                    n_triangle=8
                 elseif(check_string(buffer,"C3")) then
                    axes="C3"
                    axes_c3=.true.
                    axes_c4=.false.
                    n_triangle=8
                 elseif(check_string(buffer,"C4")) then
                    axes="C4"
                    axes_c4=.true.
                    axes_c3=.false.
                    n_triangle=8
                 elseif(check_string(buffer,"DIRECTION")) then
                    read(buffer,*,err=100) buffer1, direction(:)
                 elseif(check_string(buffer,"SCREEP_POINTS")) then
                    read(buffer,*,err=100) buffer1, n_screep_points
                 end if
              end do
           end if

10         direction(:)=direction(:) / &
                sqrt(dot_product(direction,direction))
           n_div=sqrt(float(n_screep_points/n_triangle))-1
           n_screep_points=n_triangle*(n_div+1)**2

           format="(10X, 'Screep Center : ', 3(f17.12,3x) )"
           write(output_unit,format) Rc_screep(:)
           Rc_screepAU=Rc_screep/0.529177_r8_kind
           format="(/10X, 'R_screep  : ',f10.5)"
           write(output_unit,format) r_screep
           r_screep=r_screep/0.529177_r8_kind
           format="(/10X, 'Axes type : 'a10)"
           write(output_unit,format) axes
           format="(/10X,'Orientation(x,y,z) :',3(f10.5,3x))"
           write(output_unit,format) direction(:)
           format="(/10X, 'Num.screep points :', i6)"
           write(output_unit,format) n_screep_points
           format="( 10X, 'N_div.            :', i6)"
           write(output_unit,format) n_div
           format="(/10X, 'SCREEP-procedure : ',a)"
           write(output_unit,format) trim(screep_procedure)
           write(output_unit,"(10X,'Spherical Integration scheme = ',i2)") i_scheme

           do i_atom=1,n_atoms_in_cell
              atom(i_atom)%rc=atom(i_atom)%rc-Rc_screepAU
              atom(i_atom)%rl=atom(i_atom)%rl-matmul(Rc_screepAU,b)/(2.0_r8_kind*pi)
           end do

           write(output_unit,*) "============================================================"
           format="(/34X,'X',17x,'Y',20x,'Z',17x,'Q' )"
           write(output_unit,format)
           do i_atom=1,n_atoms_in_cell
              format="(7X,i4,4x,a,3x,3(f17.12,3x),3x, f7.3 )"
              write(output_unit,format) i_atom, atom(i_atom)%name, atom(i_atom)%rc,&
                   atom(i_atom)%q
           end do
           format="(/34X,'a_1',15x,'a_2',17x,'a_3',15x,'Q')"
           write(output_unit,format)
           do i_atom=1,n_atoms_in_cell
              format="(7X,i4,4x,a,3x,3(f17.12,3x),3x,f7.3 )"
              write(output_unit,format) i_atom, atom(i_atom)%name, atom(i_atom)%rl,&
                   atom(i_atom)%q
           end do
           write(output_unit,*) "============================================================"

           return

100        call wrong_line(epeinp)

         end subroutine read_screep_sphere
         !==============================================================
         subroutine read_pc_array

           character*80 :: buffer,buffer1

           read_block=find_block(input_unit,"&PC_ARRAY",epeinp)
           if(.not.read_block) then
!              call error_handler("Read_input: No screep data found")
              goto 10
           else
              do 
                 call read_block_line(input_unit,buffer,end_block,comment_line,empty_line,epeinp)
                 if(empty_line) cycle
                 if(comment_line) cycle
                 if(end_block) exit

                 call upcase(buffer)
                 if(check_string(buffer,"RADIUS")) then
                    read(buffer,*,err=100) buffer1, R_max_pc
                 end if
              end do
           end if

10         R_max_pc = R_max_pc/0.529177_r8_kind
           format="( 10X/' PC-array size : ', f15.10, ' (a.u.)' )"
           write(output_unit,format) R_max_pc

           return

100        call wrong_line(epeinp)

         end subroutine read_pc_array
         !==============================================================
         subroutine read_ewald_parameter

           character*80 :: buffer,buffer1

           read_block=find_block(input_unit,"&EWALD_PARAMETER",epeinp)
           if(.not.read_block) then
              goto 10
           else
              do 
                 call read_block_line(input_unit,buffer,end_block,comment_line,empty_line,epeinp)
                 if(empty_line) cycle
                 if(comment_line) cycle
                 if(end_block) exit

                 call upcase(buffer)
                 if(check_string(buffer,"ETA")) then
                    read(buffer,*,err=100) buffer1, eta
                 end if
              end do
           end if

10         continue
           format="( 10X/' Ewald parameter : ', f15.10 )"
           write(output_unit,format) eta

           return

100        call wrong_line(epeinp)

         end subroutine read_ewald_parameter
         !==============================================================
         subroutine read_cluster

           integer(i4_kind), parameter :: gx_max=200_r8_kind
           character*256 :: buffer
           integer(i4_kind) :: nn,i,i_atom
           real(r8_kind) :: anm
           character*2 :: name(99)= &
                (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne", &
                "Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca", &
                "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
                "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr", &
                "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn", &
                "Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd", &
                "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tu","Yb", &
                "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg", &
                "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", &
                "Pa","U ","Np","Pu","Am","Cm","Bk","Cf","XX"/)

           allocate(cluster(gx_max), stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("Read_input: Allocation failed [5]")

           if(gx_file_exist) then
              n_cluster_atoms=0
              do
                 read(gx_unit,'(a256)',err=100) buffer
                 read(buffer,*,err=100) anm
                 if(anm <= 0.5_r8_kind) exit
                 read(buffer,*,err=100) anm,tmp_r,nn
                 if(nn == 0) cycle
                 n_cluster_atoms=n_cluster_atoms+1
                 if(n_cluster_atoms > gx_max) &
		      call error_handler("Read_input: GXFILE contains too many atoms (> 200)")
                 i=int(anm)
                 cluster(n_cluster_atoms)%name=name(i)
                 cluster(n_cluster_atoms)%rc=tmp_r*0.529177_r8_kind
                 cluster(n_cluster_atoms)%rl=matmul(tmp_r,b)*0.529177_r8_kind/(2_r8_kind*pi)
              end do
           else
              goto 200
100           write(*,*) "WARNING: GXFILE has been read in with error - no cluster used"
200           n_cluster_atoms=1
              cluster(1)%name="XX"
              cluster(1)%rc=Rc_screep
              cluster(1)%rl=matmul(Rc_screep,b)/(2_r8_kind*pi)
           end if

           deallocate(tmp_r, stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("Read_input: Deallocation failed [6]")

           write(output_unit,*) "Cluster atoms in angstrom"
           format="(/34X,'X',17x,'Y',20x,'Z')"
           write(output_unit,format)
           do i_atom=1,n_cluster_atoms
              format="(7X,i4,4x,a,3x,3(f17.12,3x) )"
              write(output_unit,format) i_atom,cluster(i_atom)%name,cluster(i_atom)%rc
              cluster(i_atom)%rc=(cluster(i_atom)%rc-Rc_screep)/0.529177_r8_kind
              cluster(i_atom)%rl=(cluster(i_atom)%rl-matmul(Rc_screep,b)/(2_r8_kind*pi))/0.529177_r8_kind
           end do

           do i_atom=1,n_cluster_atoms
              if(sqrt(dot_product((cluster(i_atom)%rc), &
                   (cluster(i_atom)%rc)) )>=r_screep) then
                 call error_handler("Error : Cluster atom is out of sphere")
              end if
           end do

           write(output_unit,*) "Cluster atoms in a.u."
           format="(/34X,'X',17x,'Y',20x,'Z')"
           write(output_unit,format)
           do i_atom=1,n_cluster_atoms
              format="(7X,i4,4x,a,3x,3(f17.12,3x) )"
              write(output_unit,format) i_atom,cluster(i_atom)%name,cluster(i_atom)%rc
           end do
           format="(/34X,'a_1',15x,'a_2',17x,'a_3')"
           write(output_unit,format)
           do i_atom=1,n_cluster_atoms
              format="(7X,i4,4x,a,3x,3(f17.12,3x) )"
              write(output_unit,format) i_atom,cluster(i_atom)%name,cluster(i_atom)%rl
           end do

         end subroutine read_cluster
         !==============================================================
         function find_block(file_id,block_name,file_num)

           logical :: find_block
           integer(i4_kind), intent(in) :: file_id,file_num
           character(len=*), intent(in) :: block_name

           character*80 :: string, message
           character*4 :: number

           rewind file_id
           line_number=0
           do
              line_number=line_number+1
              read(file_id,'(a80)', end=100, err=200) string
              if(check_comment(string)) cycle
              call upcase(string)
              if (check_string(string,block_name)) then
                 current_block=block_name
                 find_block= .true.
                 exit
              end if
           end do
           return

100        find_block= .false.
           return

200        call wrong_line(file_num)

         end function find_block
         !==============================================================
         subroutine read_block_line(file_id,line,end,comment,empty,file_num)

           integer(i4_kind), intent(in) :: file_id,file_num
           character(len=*), intent(out) :: line
           logical, intent(out) :: end,comment,empty

           character*80 :: string,message

           end=.false.; comment=.false.;empty=.false.
           line_number=line_number+1
           read(file_id,'(a80)', end=100, err=200) line
           string=line
           call upcase(string)
           if(check_empty_line(string)) empty=.true.
           if(check_comment(string)) comment=.true.
           if(check_string(string,"&END")) end=.true.
    
           return

100        message = trim("Read_"//files(file_num)//" : Unterminated block "//trim(current_block))
           call error_handler(message)

200        call wrong_line(file_num)

         end subroutine read_block_line
         !==============================================================
         function check_comment(string)

           logical :: check_comment
           character(len=*), intent(in) :: string

           check_comment=.false.
           if(index(string,"#") == 1) check_comment=.true.

         end function check_comment
         !==============================================================
         function check_empty_line(string)
           logical :: check_empty_line
           character(len=*), intent(in) :: string

           check_empty_line=.false.
           if(len(trim(string)) == 0) check_empty_line=.true.

         end function check_empty_line
         !==============================================================
         subroutine wrong_line(file_num)

           integer(i4_kind) :: file_num
           character*80 :: message
           character*4 :: number

           write(number,'(i4)') line_number
           message = trim("Read_"//files(file_num)//" : The error reading line "//trim(number))
           call error_handler(message)

         end subroutine wrong_line
         !==============================================================
         subroutine upcase(string)

           character(len=*) :: string
           integer(kind=i4_kind) :: i,ln,ich

           ln = len(trim(string))
           do i=1,ln
              ich=iachar(string(i:i))
              if(ich>=97 .and. ich<=122) then
                 ich=ich-32
                 string(i:i)=achar(ich)
              end if
           end do

         end subroutine upcase
         !==============================================================
         function check_string(string,word)

           character(len=*), intent(in) :: string
           character(len=*), intent(in) :: word

           logical :: check_string

           if(index(string,word) /= 0) then
              check_string = .true.
           else
              check_string = .false.
           end if

         end function check_string

       end subroutine ewald_input_new
!#####################################################

       subroutine screep_generate
         implicit none
         integer(kind=i4_kind), allocatable, dimension(:,:) :: index
         integer(kind=i4_kind)                              :: alloc_stat,i_point
         integer(kind=i4_kind)                              :: N_vertex,i_triangle
         integer(kind=i4_kind)                              :: i_div,step_index,screep_unit
         integer(kind=i4_kind)                              :: n_step, i_vertex,nnn
         real(kind=r8_kind), allocatable, dimension(:)      :: v_1, v_2
         real(kind=r8_kind), allocatable, dimension(:)      :: p_0,p_1,p_2,p_3
         real(kind=r8_kind), allocatable, dimension(:)      :: tmp_axes, tmp_vector
         real(kind=r8_kind), allocatable, dimension(:)      :: rot_axes 
         real(kind=r8_kind), allocatable, dimension(:,:)    :: u
         real(kind=r8_kind), allocatable, dimension(:,:)    :: mv1, mv2, mv3
         real(kind=r8_kind), dimension(3)                   :: pp_1,pp_2,pp_3
         real(kind=r8_kind)                                 :: scalar,delta,square,s_sum
         real(kind=r8_kind)                                 :: sp_square, sp_sum
         real(kind=r8_kind)                                 :: mean, mean2
         character(len=100)                                 :: format,screep_file

         screep_file=trim(work_dir)//"/"//"screep.xyz"     
         screep_unit=8      ! to be changed by "open_get ..."
         open(screep_unit,file=trim(screep_file), status="unknown")

         allocate(screep_point(n_screep_points), stat=alloc_stat)
         if(alloc_stat.ne.0) &
              call error_handler("Screep_generate [1]: Allocation failed")
       
         allocate(index(n_triangle,3),stat=alloc_stat)
         if(alloc_stat.ne.0) &
              call error_handler("Screep_generate [2]: Allocation failed")
         index=reshape((/1,2,3,4,1,2,3,4,2,3,4,1, &
                            2,3,4,1,5,5,5,5,6,6,6,6/),(/n_triangle,3/))
         N_vertex=(maxval(index))
         allocate(vertex(3,N_vertex), stat=alloc_stat)
         if(alloc_stat.ne.0) &
              call error_handler("Screep_generate [3]: Allocation failed")
!         C_4 || z-axes (default, will be changed depending on axes chosen
            vertex=reshape((/ 1.0_r8_kind, 0.0_r8_kind, 0.0_r8_kind,    &
                              0.0_r8_kind, 1.0_r8_kind, 0.0_r8_kind,    &
                             -1.0_r8_kind, 0.0_r8_kind, 0.0_r8_kind,    &
                              0.0_r8_kind,-1.0_r8_kind, 0.0_r8_kind,    &
                              0.0_r8_kind, 0.0_r8_kind, 1.0_r8_kind,    &
                              0.0_r8_kind, 0.0_r8_kind,-1.0_r8_kind /), &
                             (/3,N_vertex/) )
          if(.not.axes_c3 .and. .not.axes_c4 ) &
             call error_handler("Screep_generate [4] : Axes type not defined ")
          allocate(tmp_axes(3), tmp_vector(3), stat=alloc_stat)
          if(alloc_stat.ne.0) &
               call error_handler("Screep_generate [5]: Allocation failed")
          if(axes_c4) tmp_axes=(/0.0_r8_kind, 0.0_r8_kind, 1.0_r8_kind/)
          if(axes_c3) then
             tmp_axes(:)=(vertex(:,1)+vertex(:,2)+vertex(:,5))/3.0_r8_kind
             scalar=sqrt( dot_product(tmp_axes,tmp_axes) )
             tmp_axes(:)=tmp_axes(:)/scalar
          end if

          allocate(rot_axes(3), u(3,3), stat=alloc_stat)
          if(alloc_stat.ne.0) &
               call error_handler("Screep_generate [6]: Allocation failed")
          delta=acos( dot_product(direction,tmp_axes) / &
                     (sqrt(dot_product(direction,direction)) * &
                      sqrt(dot_product(tmp_axes,tmp_axes)) ) )
          if(abs(delta)/=0.0_r8_kind) then
             rot_axes=cross_product(tmp_axes,direction)
             rot_axes=rot_axes/sqrt(dot_product(rot_axes,rot_axes))
             call rot(u,rot_axes,delta)
             do i_vertex=1,n_vertex
                tmp_vector=matmul(u,vertex(:,i_vertex))
                vertex(:,i_vertex)=tmp_vector(:)
             end do
          end if
	if(symmetry_c3v) then
	   tmp_axes=(/2.9511013 , -0.7907452, 0.0/)
	   rot_axes=(/0.0, 1.0, 0.0/)
	   delta=acos( dot_product(tmp_axes, rot_axes)/&
		       sqrt(dot_product(tmp_axes,tmp_axes)))
	   rot_axes=(/0.0,0.0,1.0/)
	   call rot(u,rot_axes,delta)
	   do i_vertex=1,n_vertex
		tmp_vector=matmul(u,vertex(:,i_vertex))
		vertex(:,i_vertex)=tmp_vector(:)
           end do
	endif

          deallocate(tmp_vector,tmp_axes,rot_axes, stat=alloc_stat)
          if(alloc_stat.ne.0) &
               call error_handler("Screep_generate [7]: deallocation failed")
!
!                          Fill the triangle-faces
!
        allocate(v_1(3),v_2(3),p_0(3),p_1(3),p_2(3),p_3(3), stat=alloc_stat)
        if(alloc_stat.ne.0) &
             call error_handler("Screep_generate [8]: Allocation failed")
        allocate(mv1(3,n_screep_points),mv2(3,n_screep_points),mv3(3,n_screep_points))
        i_point=0
        triangle:  do i_triangle=1,n_triangle
           v_1(:)=(vertex(:,index(i_triangle,2))-vertex(:,index(i_triangle,1)))/&
                  (n_div+1)
           v_2(:)=(vertex(:,index(i_triangle,3))-vertex(:,index(i_triangle,1)))/&
                  (n_div+1)

        do i_div=1,n_div+1
           p_0(:)=vertex(:,index(i_triangle,1))+add_vector(0,v_1,i_div-1,v_2,3)
           step_index=-1.
           do n_step=1,2*(n_div-i_div)+3
              select case(step_index)
                 case(-1)
                    p_1=p_0
                    p_2=p_0+add_vector(1,v_1, 0,v_2,3)
                    p_3=p_0+add_vector(0,v_1, 1,v_2,3)
                    p_0=p_3
                    i_point=i_point+1
                    screep_point(i_point)%rc=(p_1+p_2+p_3) / 3.0_r8_kind
                    step_index=-step_index
                    mv1(:,i_point)=p_1
                    mv2(:,i_point)=p_2
                    mv3(:,i_point)=p_3
                 case(1)
                    p_1=p_0
                    p_2=p_0+add_vector(1,v_1, 0,v_2,3)
                    p_3=p_0+add_vector(1,v_1,-1,v_2,3)
                    p_0=p_3
                    i_point=i_point+1
                    screep_point(i_point)%rc=(p_1+p_2+p_3) / 3.0_r8_kind
                    step_index=-step_index
                    mv1(:,i_point)=p_1(:)
                    mv2(:,i_point)=p_2(:)
                    mv3(:,i_point)=p_3(:)
              end select
           end do
        end do

        end do triangle
        deallocate(v_1,v_2,p_0,p_1,p_2,p_3, stat=alloc_stat)
        if(alloc_stat.ne.0) &
             call error_handler("Screep_generate [9]: Deallocation failed")

        s_sum=0.0_r8_kind
        mean=0.0_r8_kind
        mean2=0.0_r8_kind
        do i_point=1,n_screep_points
           scalar=sqrt(dot_product(screep_point(i_point)%rc,screep_point(i_point)%rc))
           screep_point(i_point)%rc=R_screep*screep_point(i_point)%rc / scalar
!-----------------
           mv1(:,i_point)=R_screep*mv1(:,i_point) / &
                          sqrt(dot_product(mv1(:,i_point),mv1(:,i_point)))
           mv2(:,i_point)=R_screep*mv2(:,i_point) / &
                          sqrt(dot_product(mv2(:,i_point),mv2(:,i_point)))
           mv3(:,i_point)=R_screep*mv3(:,i_point) / &
                          sqrt(dot_product(mv3(:,i_point),mv3(:,i_point)))
           square =sqrt(abs(dot_product(cross_product(mv3(:,i_point)-mv1(:,i_point),  &
                                                      mv3(:,i_point)-mv2(:,i_point)), &
                                        cross_product(mv3(:,i_point)-mv1(:,i_point),  &
                                                      mv3(:,i_point)-mv2(:,i_point)))))/2

           sp_square=sph_square(mv1(:,i_point),mv2(:,i_point),mv3(:,i_point),R_screep)
           screep_point(i_point)%square=sp_square
           s_sum  = s_sum + square
           sp_sum = sp_sum + sp_square
           mean=mean + sp_square
           mean2=mean2 + sp_square**2
           format="(10X/i4,5x,'S=',f10.5,5x,'Sph=',f10.5)"

           write(output_unit,format) i_point,square,sp_square
!-----------------
        end do
        write(output_unit,*) s_sum, sp_sum, 4*3.14159265358_r8_kind*R_screep**2
        write(output_unit,*) "Dispersion : ", &
                           sqrt(mean2/n_screep_points - (mean/n_screep_points)**2)
        deallocate(mv1,mv2,mv3,vertex, stat=alloc_stat)
        if(alloc_stat.ne.0) &
             call error_handler("Screep_generate [10]: Deallocation failed")
!                    Save generated screep_points in MOLDEN ".xyz" format 
        write(screep_unit,*) n_screep_points
        write(screep_unit,*) "comment line"
        do i_point=1,n_screep_points
           format="('H',3(5x,f12.7))"
         write(screep_unit,format) &
	       0.529177_r8_kind*(screep_point(i_point)%rc(:)+Rc_screepAU)
        end do
        close(screep_unit, status='keep')

       end subroutine screep_generate

       function add_vector(n_1,vec_1,n_2,vec_2,dim)
         use type_module
         implicit none
         integer(kind=i4_kind), intent(in)            :: n_1,n_2,dim
         real(kind=r8_kind), dimension(3), intent(in) :: vec_1,vec_2
         real(kind=r8_kind), dimension(dim)             :: add_vector

         add_vector(:)=n_1*vec_1(:) + n_2*vec_2(:)
         return
       end function add_vector

       function cross_product(v_1,v_2)
         real(kind=8)               :: cross_product(3)
         real(kind=8), dimension(3) :: v_1, v_2

         cross_product(1)=v_1(2)*v_2(3)-v_1(3)*v_2(2)
         cross_product(2)=v_1(3)*v_2(1)-v_1(1)*v_2(3)
         cross_product(3)=v_1(1)*v_2(2)-v_1(2)*v_2(1)
         return
       end function cross_product

!======================================================================!
!                  rotation around axes (nx,ny,nz) on d                !
!======================================================================!
      subroutine rot(u,rot_axes,d)
        implicit none
        real(kind=r8_kind),intent(inout), dimension(3)   :: rot_axes(3)
        real(kind=r8_kind),intent(inout), dimension(3,3) :: u
        real(kind=r8_kind),intent(in)                    :: d
        real(kind=r8_kind)                               :: cosd,sind,nx,ny,nz

        nx=rot_axes(1)
        ny=rot_axes(2)
        nz=rot_axes(3)
        cosd=cos(d)
        sind=sin(d)
        u(1,1)=cosd + (1.d0-cosd)*nx*nx
        u(2,1)=(1.d0-cosd)*nx*ny + sind*nz
        u(3,1)=(1.d0-cosd)*nx*nz - sind*ny
        u(1,2)=(1.d0-cosd)*nx*ny -sind*nz
        u(2,2)=cosd + (1.d0-cosd)*ny*ny
        u(3,2)=(1.d0-cosd)*ny*nz + sind*nx
        u(1,3)=(1.d0-cosd)*nx*nz + sind*ny
        u(2,3)=(1.d0-cosd)*ny*nz - sind*nx
        u(3,3)=cosd + (1.d0-cosd)*nz*nz
        
      end subroutine rot

     function sph_square(v_1,v_2,v_3,R)
        implicit none
        real(kind=r8_kind)                           :: sph_square
        real(kind=r8_kind), intent(in), dimension(3) :: v_1,v_2,v_3
        real(kind=r8_kind), intent(in)               :: R
        real(kind=r8_kind)                           :: alpha, beta, gamma, pi

        pi=acos(-1.0_r8_kind)
        alpha=acos( dot_product(cross_product(v_1,v_3),cross_product(v_2,v_3))/ &
                    (sqrt(dot_product(cross_product(v_1,v_3),cross_product(v_1,v_3))) * &
                     sqrt(dot_product(cross_product(v_2,v_3),cross_product(v_2,v_3))) ) )
        beta =acos( dot_product(cross_product(v_3,v_1),cross_product(v_2,v_1))/ &
                    (sqrt(dot_product(cross_product(v_3,v_1),cross_product(v_3,v_1))) * &
                     sqrt(dot_product(cross_product(v_2,v_1),cross_product(v_2,v_1))) ) )
        gamma=acos( dot_product(cross_product(v_3,v_2),cross_product(v_1,v_2))/ &
                    (sqrt(dot_product(cross_product(v_3,v_2),cross_product(v_3,v_2))) * &
                     sqrt(dot_product(cross_product(v_1,v_2),cross_product(v_1,v_2))) ) )
        sph_square= R**2 *(alpha+beta+gamma-pi)

      end function sph_square

!===================================================================================!
!                                                                                   !
!       Purpose: subroutine generates the set of lattice vectors,                   !
!                both in real and reciprocal space;                                 !
!                In case of rec.space, the variables mu1, mu2, mu3                  !
!                define lattice coordinates of k-vector (or k-point)                !
!                If option "full_print" is included the full output will            !
!                be generated (i.e. vector, cartesian coordiates and                !
!                modules of lattice vetors)                                         !
!
      subroutine generate_lattice(space,mu_1,mu_2,mu_3,r_max,print_option,full_print)
        implicit none
        real(kind=8),      intent(in)      :: mu_1, mu_2, mu_3
        real(kind=8),      intent(in)      :: r_max
        logical,           intent(in)      :: print_option
        logical, optional, intent(in)      :: full_print
        character(len=*),  intent(in)      :: space ! 'direct' or 'reciprocal'

!                              local variables

        real(kind=r8_kind), allocatable, dimension(:)      :: aux_vector,R,h
        real(kind=r8_kind), allocatable, dimension(:)      :: v_temp_cartesian, v_temp_vector
        real(kind=r8_kind), allocatable, dimension(:,:)    :: a_lat,v
        real(kind=r8_kind)                                 :: a_min, mod_temp 
        integer(kind=i4_kind)                              :: i,k,l
        integer(kind=i4_kind)                              :: num_slab, num_in_slab
        integer(kind=i4_kind)                              :: num_omitted,n_1, n_2, n_3
        integer(kind=i4_kind)                              :: nu, i_cycle, num, m, j
        integer(kind=i4_kind)                              :: i_vector, j_vector,i_print
        integer(kind=i4_kind)                              :: i_first, i_second, alloc_stat
        logical                                            :: in_slab, empty_cycle
        character(len=32) print_format

        allocate(a_lat(3,3))
        select case(space)
        case('direct')
           a_lat = a
        case('reciprocal')
           a_lat = b
        end select

!      First, empty cycle to obtain the number of vectors (num_vector)
!                          Inscribe a "R_max sphere"
        allocate(R(3),v(3,3),h(6), stat=alloc_stat)
        if(alloc_stat.ne.0) &
             call error_handler("generate_lattice [1]: Allocation failed")
        num_slab=1
        h=0.0_r8_kind
        do while(.not.all(h>r_max))
           do l=1,2
              do i=1,3
                 j=i+1
                 if(j.gt.3) j=j-3
                 k=j+1
                 if(k.gt.3) k=k-3
                 if(l==1) then
                    v(:,j)=2*num_slab*a_lat(:,j)
                    v(:,k)=2*num_slab*a_lat(:,k)
                    R(:)= (mu_1-num_slab)*a_lat(:,1) + &
                          (mu_2-num_slab)*a_lat(:,2) + (mu_3-num_slab)*a_lat(:,3)
                 else
                    v(:,j)=-2*num_slab*a_lat(:,j)
                    v(:,k)=-2*num_slab*a_lat(:,k)
                    R(:)=-(mu_1+num_slab)*a_lat(:,1) - &
                          (mu_2+num_slab)*a_lat(:,2) - (mu_3+num_slab)*a_lat(:,3)
                 end if
                 h((l-1)*3+i)= abs(dot_product(R,cross_product(v(:,j),v(:,k)))) / &
                               sqrt(dot_product(cross_product(v(:,j),v(:,k)),&
                                    cross_product(v(:,j),v(:,k))))
              end do
           end do
           num_slab=num_slab+1
        end do
        deallocate(R,v,h, stat=alloc_stat)
        if(alloc_stat.ne.0) &
             call error_handler("generate_lattice [1]: Deallocation failed")

        allocate(aux_vector(3) )
        loop_empty: do i_cycle=1,2

        empty_cycle=i_cycle==1 
        if(empty_cycle) then
           if(associated(lattice%v_cartesian)) deallocate(lattice%v_cartesian)
           if(associated(lattice%v_vector))    deallocate(lattice%v_vector)
           if(associated(lattice%v_mod))       deallocate(lattice%v_mod)      
        else
           allocate( lattice%v_cartesian(lattice%num_vector,3), & 
                 lattice%v_vector(lattice%num_vector,3),    & 
                 lattice%v_mod(lattice%num_vector) )
        end if
        lattice%num_vector=0
        do n_1 = -num_slab, num_slab
           do n_2 = -num_slab, num_slab
              do n_3 = -num_slab, num_slab
                    num_in_slab=num_in_slab+1
                    aux_vector(:) = (n_1 + mu_1) * a_lat(:,1)  +  &
                                    (n_2 + mu_2) * a_lat(:,2)  +  &
                                    (n_3 + mu_3) * a_lat(:,3)
                    if( sqrt(dot_product(aux_vector,aux_vector)) <= r_max) then
                       lattice%num_vector=lattice%num_vector+1
                       if(.not.empty_cycle) then
                          lattice%v_cartesian(lattice%num_vector,:) = aux_vector(:)
                          lattice%v_vector(lattice%num_vector,1) = n_1
                          lattice%v_vector(lattice%num_vector,2) = n_2
                          lattice%v_vector(lattice%num_vector,3) = n_3
                          lattice%v_mod(lattice%num_vector)      = sqrt( &
                               dot_product(aux_vector,aux_vector) )
                       end if
                    end if
                 end do
          end do
       end do
  
 end do loop_empty

        deallocate( aux_vector, a_lat )
        allocate(v_temp_cartesian(3), v_temp_vector(3) )
        do i_vector=1, lattice%num_vector
           a_min = lattice%v_mod(i_vector)
           num=i_vector

           do j_vector=i_vector, lattice%num_vector
              if(lattice%v_mod(j_vector) <= a_min ) then
                 a_min = lattice%v_mod(j_vector)
                 num=j_vector
              end if
           end do ! j_vector

           v_temp_cartesian(:) = lattice%v_cartesian(num,:)
           v_temp_vector(:)    = lattice%v_vector(num,:)  
           mod_temp            = lattice%v_mod(num)
           lattice%v_cartesian(num,:)      = lattice%v_cartesian(i_vector,:)
           lattice%v_cartesian(i_vector,:) = v_temp_cartesian(:)
           lattice%v_vector(num,:)         = lattice%v_vector(i_vector,:)
           lattice%v_vector(i_vector,:)    = v_temp_vector(:)
           lattice%v_mod(num)              = lattice%v_mod(i_vector)
           lattice%v_mod(i_vector)         = mod_temp

        end do ! i_vector
        deallocate(v_temp_cartesian, v_temp_vector )


        do i_cycle=1,2
           empty_cycle=i_cycle==1
           lattice%n_spheres=1
           if(.not.empty_cycle) lattice%num_in_sphere(lattice%n_spheres)=1 
           do i_vector=2,lattice%num_vector
              if( abs(lattice%v_mod(i_vector)-lattice%v_mod(i_vector-1)) <= 1.d-3 ) then
                 if(.not.empty_cycle) lattice%num_in_sphere(lattice%n_spheres)=&
                                      lattice%num_in_sphere(lattice%n_spheres)+1
              else 
                 lattice%n_spheres=lattice%n_spheres + 1
                 if(.not.empty_cycle) lattice%num_in_sphere(lattice%n_spheres)=1
              end if
           end do
           if(empty_cycle) allocate(lattice%num_in_sphere(lattice%n_spheres) )
        end do

        print_format="(10x, 'Number of vectors : ', i5/10x, 'Number of shells :',i5 )"
        write(output_unit,*) lattice%num_vector, lattice%n_spheres
        do i_print=1,(lattice%n_spheres-1)/10+1
           print_format="(10x,10i5)"
           write(output_unit,print_format) &
                (lattice%num_in_sphere(m),m=((i_print-1)*10+1),&
                                              min(i_print*10,lattice%n_spheres))
        end do
              
        if(print_option) then
           if( present(full_print) ) then
              write(6,'(/20x,"Coordinates in vector units"/)')
              do i_print=1,lattice%num_vector
                 i_first=(i_print-1)*5 + 1
                 i_second=i_print*5
                 if(i_second > lattice%num_vector) i_second=lattice%num_vector
                 print_format='(3x, 5("(",3i4,")",1x), i3)'
                 write(6,print_format) ( (lattice%v_vector(m,j),j=1,3), &
                                      m=i_first,i_second),i_second 
                 if(i_second==lattice%num_vector) exit               
              end do

        write(6,'(/20x,"Cartesian coordinates"/)')
        do i_print=1,lattice%num_vector
           i_first=(i_print-1)*5 + 1
           i_second=i_print*5
           if(i_second > lattice%num_vector) i_second=lattice%num_vector
           print_format='( 3x, 5("(",3f9.4,")",2x), / )'
           write(6,print_format) ( (lattice%v_cartesian(m,j),j=1,3), &
                                                 m=i_first,i_second) 
           if(i_second==lattice%num_vector) exit 
        end do
     end if

        write(6,'(/20x, "Vector modules"/)')
        do i_print=1,lattice%num_vector
           i_first=(i_print-1)*5 + 1
           i_second=i_print*5
           if(i_second > lattice%num_vector) i_second=lattice%num_vector
           print_format='(3x, 5g12.5, i5)'
           write(6,print_format) (lattice%v_mod(m), m=i_first,i_second), i_second
           if(i_second==lattice%num_vector) exit 
        end do
     end if ! print_option
   end subroutine generate_lattice


      subroutine pc_generate
        implicit none
        type(atom_type),    pointer, dimension(:)     :: atom_wz, tmp
        type(point_charge), pointer, dimension(:)     :: tmp_pc
        real(kind=r8_kind), allocatable, dimension(:) :: rs
        real(kind=r8_kind)                            :: mu_1,mu_2,mu_3,q_cell,q_tot,r_max,R
        real(kind=r8_kind),              dimension(3) :: p_cell
        integer(kind=i4_kind)                         :: N_wz,N_fragment_atoms,i_atom,n_equiv,i,j,k
        integer(kind=i4_kind)                         :: alloc_stat, scell_unit,pc_xyz_unit
        character(len=100)                            :: format,scell_file,pc_xyz_file

        scell_file=trim(work_dir)//"/"//"cell.xyz"
        scell_unit=7                                    ! to be changed by "openget ..."
        open(scell_unit, file=trim(scell_file), status="unknown")
        pc_xyz_file=trim(work_dir)//"/"//"Ewald_PC.xyz"
        pc_xyz_unit=10                                  ! to be changed by "openget ..."
        open(pc_xyz_unit, file=trim(pc_xyz_file), status="unknown")
!                                  Form symmetric cell
        N_wz=0
        do i_atom=1,n_atoms_in_cell
           mu_1=atom(i_atom)%rl(1)
           mu_2=atom(i_atom)%rl(2)
           mu_3=atom(i_atom)%rl(3)
           r_max=sqrt( dot_product(atom(i_atom)%rc,atom(i_atom)%rc) )
           if(r_max==0.0_r8_kind) r_max=1.0_r8_kind
           do 
              call generate_lattice("direct",mu_1,mu_2,mu_3,r_max,print_vectors)
              if(lattice%n_spheres >= 2) exit
              r_max=2*r_max
           end do        !          First shell only
           write(output_unit,*) 'atom ',i_atom, 'n_equiv= ',lattice%num_in_sphere(1) 
           n_equiv=lattice%num_in_sphere(1)
           if(i_atom > 1) then
              allocate(tmp(N_wz), stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("pc_generate [1]: Allocation failed")
              tmp=atom_wz
              deallocate(atom_wz, stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("pc_generate [2]: Deallocation failed")
              allocate(atom_wz(N_wz+n_equiv), stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("pc_generate [3]: Allocation failed")
              do i=1,N_wz
                 atom_wz(i)=tmp(i)
              end do
              deallocate(tmp, stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("pc_generate [4]: Deallocation failed")
           else
              allocate(atom_wz(n_equiv), stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("pc_generate [5]: Allocation failed")
           end if
           do i=1,n_equiv
              atom_wz(N_wz+i)%rc(:)=lattice%v_cartesian(i,:)
              atom_wz(N_wz+i)%rl(:)=atom(i_atom)%rl(:)+lattice%v_vector(i,:)
              atom_wz(N_wz+i)%q=atom(i_atom)%q / n_equiv
              atom_wz(N_wz+i)%name=atom(i_atom)%name
           end do
           N_wz=N_wz+n_equiv
        end do ! i_atom
!                                Control cell-charge and store sym.cell
        q_cell=0.0_r8_kind
        p_cell=0.0_r8_kind
        write(output_unit,"(/37(1h#),' Symmetric cell ',37(1h#)  )" )
        write(output_unit,"(10X,'Symmetric cell is stored in : ',a)") scell_file
        write(scell_unit,"(i7/'comment line')" ) N_wz        
        write(output_unit,"(/34X,'X',17x,'Y',20x,'Z',17x,'Q' )")
        do i_atom=1,N_wz
           q_cell=q_cell + atom_wz(i_atom)%q
           p_cell(:)=p_cell(:) + atom_wz(i_atom)%rc(:)*atom_wz(i_atom)%q
           format="(7X,i4,4x,a,3x,3(f17.12,3x),3x, f7.3 )" 
           write(output_unit,format) i_atom, atom_wz(i_atom)%name, atom_wz(i_atom)%rc(:),&
                                             atom_wz(i_atom)%q
           format="(a4,4x,3(f12.7,3x),5x,f10.7,3x,i4)"
           write(scell_unit,format) atom_wz(i_atom)%name,atom_wz(i_atom)%rc(:), & 
                                    atom_wz(i_atom)%q, i_atom
        end do
        format="(10x, 'Charge of symmetric cell :', f10.5 &
               &/10x, 'Dipole of symmetric cell :', 3(f10.5,3x)/90(1h#))"
        write(output_unit,format) q_cell,p_cell(:)
!                             Make PC_array around the sym.cell
        call generate_lattice("direct",0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,R_max_pc,&
                               print_vectors,full_print)
        allocate(rs(3), stat=alloc_stat)
        if(alloc_stat.ne.0) &
             call error_handler("pc_generate [6]: Allocation failed")
        N_ewald_pc=0
        do i_atom=1,N_wz
           lat: do i=1,lattice%num_vector
              rs(:)=atom_wz(i_atom)%rc(:)+lattice%v_cartesian(i,:)
              if(N_ewald_pc >= 1) then
                 do j=1,N_ewald_pc
!                    if( all(rs==pc(j)%rc) ) then
                     if( all( abs(rs-pc(j)%rc)<=1.d-10) ) then
                       pc(j)%q=pc(j)%q + atom_wz(i_atom)%q
                       cycle lat
                    end if
                 end do
              end if
              if(N_ewald_pc >= 1) then
                 allocate(tmp_pc(N_ewald_pc), stat=alloc_stat)
                 if(alloc_stat.ne.0) &
                      call error_handler("pc_generate [7]: Allocation failed")
                 tmp_pc=pc
                 deallocate(pc, stat=alloc_stat)
                 if(alloc_stat.ne.0) &
                      call error_handler("pc_generate [8]: Deallocation failed")
                 allocate(pc(N_ewald_pc+1),stat=alloc_stat)
                 if(alloc_stat.ne.0) &
                      call error_handler("pc_generate [9]: Allocation failed")  
                 do k=1,N_ewald_pc
                    pc(k)=tmp_pc(k)
                 end do
                 deallocate(tmp_pc, stat=alloc_stat)
                 if(alloc_stat.ne.0) &
                      call error_handler("pc_generate [10]: Deallocation failed")
              else
                 allocate(pc(1),stat=alloc_stat)
                 if(alloc_stat.ne.0) &
                      call error_handler("pc_generate [11]: Allocation failed")
              end if
              N_ewald_pc=N_ewald_pc+1
              pc(N_ewald_pc)%rc=rs
              pc(N_ewald_pc)%q=atom_wz(i_atom)%q
              pc(N_ewald_pc)%name=atom_wz(i_atom)%name
           end do lat
        end do
!                        Truncate (if needed) and store crystal fragment generated        
        if(make_crystal_fragment) then
           allocate(tmp_pc(N_ewald_pc),stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("pc_generate [12]: Allocation failed")
           N_fragment_atoms=0
           do i_atom=1,N_ewald_pc
              if( sqrt( dot_product(pc(i_atom)%rc,pc(i_atom)%rc) ) <= R_max_fragment ) then
                 N_fragment_atoms=N_fragment_atoms + 1
                 tmp_pc(N_fragment_atoms)=pc(i_atom)
              end if
           end do
           do i_atom=1,N_fragment_atoms
              pc(i_atom)=tmp_pc(i_atom)
           end do
           deallocate(tmp_pc, stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("pc_generate [13]: Deallocation failed")
           write(output_unit,"(/37(1h#),' Crystal fragment ',37(1h#)  )" )
           write(output_unit,"(10X,'Crystal fragment is stored in : ',a)") pc_xyz_file  ! repeat !!!
           write(pc_xyz_unit,"(i7/'comment line')" ) N_fragment_atoms
           write(output_unit,"(/34X,'X',17x,'Y',20x,'Z',17x,'Q' )")
           do i=1,N_fragment_atoms
              format="(7X,i4,4x,a,3x,3(f17.12,3x),3x, f7.3 )" 
              write(output_unit,format) i, pc(i)%name, pc(i)%rc(:), pc(i)%q
              format="(a4,4x,3(f17.12,3x),5x,f10.7,3x,i4)"
              write(pc_xyz_unit,format)  &
                  pc(i)%name, &
                  0.529177_r8_kind*(pc(i)%rc(:)+rc_screepAU),pc(i)%q, i  
           end do
           return
        end if ! make_crystal_fragment
!                                   Store Ewald_PC.xyz
        write(output_unit,"(/37(1h#),' Ewald PC-array ',37(1h#)  )" )
        write(output_unit,"(10X,'Ewald PC-array is stored in : ',a)") pc_xyz_file
        write(pc_xyz_unit,"(i7/'comment line')" ) N_ewald_pc
        write(output_unit,"(/34X,'X',17x,'Y',20x,'Z',17x,'Q' )")
        cluster(:)%V_pc=0.0_r8_kind
        q_tot=0.0_r8_kind
        do i=1,N_ewald_pc
           q_tot=q_tot+pc(i)%q
           format="(7X,i4,4x,a,3x,3(f17.12,3x),3x, f7.3 )" 
           write(output_unit,format)  &
                i, pc(i)%name, pc(i)%rc(:), pc(i)%q
           format="(a4,4x,3(f12.7,3x),5x,f10.7,3x,i4)"
           write(pc_xyz_unit,format) &
                pc(i)%name,0.529177_r8_kind*(pc(i)%rc(:)+rc_screepAU),pc(i)%q, i      
           do i_atom=1,n_cluster_atoms
              if( all(abs(pc(i)%rc-cluster(i_atom)%rc)<=1.d-1) ) cycle
              cluster(i_atom)%V_pc = cluster(i_atom)%V_pc + pc(i)%q / &
                   sqrt(dot_product(pc(i)%rc-cluster(i_atom)%rc,pc(i)%rc-cluster(i_atom)%rc))
           end do
        end do

        write(output_unit,"(/10x,'Total PC charge                  : ',f12.7)") q_tot
        write(output_unit,"(/10X,'Coulomb potential on cluster atoms')")
        format="(10x,i3,3x,a,4x,f12.7,4x,f12.7)"
        do i_atom=1,n_cluster_atoms
           write(output_unit,format) i_atom,cluster(i_atom)%name, cluster(i_atom)%V_pc
        end do
!!!!!!!!!!!!!!!!!!!
        do i_atom=1,n_cluster_atoms
           cluster(i_atom)%V_pc=0.0_r8_kind
           k=0
           do i=1,N_ewald_pc
              R=sqrt(dot_product(pc(i)%rc-cluster(i_atom)%rc,pc(i)%rc-cluster(i_atom)%rc))
              if( R <=1.d-1 ) cycle
              k=k+1
              write(output_unit,*) "R, i_atom, i, pc(i)%q / R  ",R, i_atom, i, pc(i)%q / R
              cluster(i_atom)%V_pc=cluster(i_atom)%V_pc + pc(i)%q / R
           end do
           write(output_unit,*) "Cl  :", i_atom, cluster(i_atom)%V_pc, N_ewald_pc,k
        end do
!!!!!!!!!!!!!!!!!!
     
      end subroutine pc_generate

      subroutine Madelung(rl_in,V_madelung,n_centers)
        implicit none
        integer(kind=i4_kind), intent(in)                       :: n_centers
        real(kind=r8_kind),              dimension(n_centers,3) :: rl_in
        real(kind=r8_kind), intent(out), dimension(n_centers)   :: V_madelung
        real(kind=r8_kind),              dimension(3)           :: rc,v
        real(kind=r8_kind)                                      :: zero=0.0_r8_kind
        real(kind=r8_kind)                                      :: x_eps,R_max,G_max
        real(kind=r8_kind)                                      :: mu_1,mu_2,mu_3
        real(kind=r8_kind)                                      :: R,G,Gr,M,pi,eps=1.d-15
        integer(kind=i4_kind)                                   :: i_atom,i_vector,i_center,i

        pi=acos(-1.d0)
        x_eps=r_f(erfc,eps,zero)
        R_max=x_eps/eta                     ! upper estimate for "direct" lattice summation
        x_eps=r_f(expf,4*eta**2*eps,0.5_r8_kind)
        G_max=2*eta*x_eps                   ! upper estimate for "direct" lattice summation
        write(output_unit,*) 'Direct Sum R_max :', R_max
        write(output_unit,*) 'Recipr.Sum G_max :', G_max, exp(-(G_max/(2*eta))**2)/G_max**2
!
        call generate_lattice("direct",zero,zero,zero,R_max,print_vectors)
        do i_center=1,n_centers
!print*, 'R', i_center
           do i=1,3
              do while(abs(rl_in(i_center,i)) > 1.0_r8_kind)
                 rl_in(i_center,i)=rl_in(i_center,i)-rl_in(i_center,i)/abs(rl_in(i_center,i))
              end do 
           end do
           rc=matmul(a,rl_in(i_center,:))
           V_madelung(i_center)=zero
           do i_atom=1,n_atoms_in_cell
              M=zero
              do i_vector=1,lattice%num_vector                 
                 v(:)=lattice%v_cartesian(i_vector,:)+atom(i_atom)%rc(:) - rc(:)  
                 R=sqrt(dot_product(v,v))
                 if(R<=1.d-1) R=zero
                 if(R/=zero) then
                    M = M + erfc(eta*R)/R
                 else
                    M = M - 2*eta/sqrt(pi)
                 end if
              end do
              V_madelung(i_center)=V_madelung(i_center) + atom(i_atom)%q * M           
           end do
        end do
!                            Add the summ over reciprocal lattice

        call generate_lattice("reciprocal",zero,zero,zero,G_max,print_vectors,full_print)
        do i_center=1,n_centers
!print*, 'G',i_center
           do i_atom=1,n_atoms_in_cell
              M=zero
              do i_vector=2,lattice%num_vector
                 G=lattice%v_mod(i_vector)
                 Gr=2*pi*dot_product(lattice%v_vector(i_vector,:), &
                                     rl_in(i_center,:)-atom(i_atom)%rl(:))
                 M = M + exp(-(G/(2*eta))**2 )*cos(Gr) / G**2
              end do
              V_madelung(i_center)=V_madelung(i_center) + &
                                   atom(i_atom)%q * (4*pi*M/V_cell - pi/(eta*eta*V_cell) )
           end do
        end do

      end subroutine Madelung

      function r_f(f,y,x0)  result(x)
        implicit none
        real(kind=8), intent(in)       :: y,x0
        real(kind=8)                   :: x,dx,x1,x2,tol=1.d-15
        real(kind=8), external         :: f

        dx=1.0_8
        x1=x0
        x2=x1+dx
        if(f(x1)>y) then
           do while(f(x2)>y)
              x1=x2
              x2=x1+dx
           end do
        else
           do while(f(x1)>y)
              x2=x1
              x1=x2-dx
           end do
        end if
        x=(x1+x2)/2.0_8
        do while( dx > tol )
           if(f(x)> y) then
              x1=x
           else
              x2=x
           end if
           x=(x1+x2)/2.0_8
           dx=abs(x2-x1)
        end do
write(output_unit,*) 'eps  =', y
write(output_unit,*) 'F(X0)=',f(x),x
      end function r_f

      function expf(x)
        real(kind=r8_kind)  :: expf,x
        expf=exp(-x**2)/x**2
      end function expf

      function erfc(x)
        !
        ! COMPLEMETARY ERROR FUNCTION USING CHEBYSHEV APPROXIMATION (YL LUKE
        ! 1975 PP123-4) D,DD,SV AS IN PRESS, NUM. REC. 1988 "CHEBEV"
        !
        ! See also modules/fermi_module.f90, epe_dir/culon_module.f90,
        ! utilities/ewald.f90
        !
        implicit none
        integer(kind=4), parameter :: na=25, nc=22
        real(kind=8), dimension(0:na) :: a
        real(kind=8), dimension(0:nc) :: c
        real(kind=8)                  :: zero=0.0, half=0.5, one=1.0,     &
                                         two=2.0,three=3.0, four=4.0,     &
                                         sqpi2=1.128379167095513
        real(kind=8)                  :: erfc,x,z,d,dd,alpha,sv
        integer(kind=4)               :: j
        DATA A /                                                          &
        .109547129977762D+1, -.289175401126989D+0,  .110456398633795D+0,  &
       -.412531882278565D-1,  .140828380706516D-1, -.432929544743143D-2,  &
        .119827190159228D-2, -.299972962353249D-3,  .683258603788747D-4,  &
       -.142469884548677D-4,  .273540877283989D-5, -.048619128719754D-5,  &
        .008038727621172D-5, -.001241841831213D-5,  .000179953258879D-5,  &
       -.000024547948775D-5,  .000003162508603D-5, -.000000385902200D-5,  &
        .000000044720291D-5, -.000000004933613D-5,  .000000000519303D-5,  &
       -.000000000052258D-5,  .000000000005037D-5, -.000000000000466D-5,  &
        .000000000000041D-5, -.000000000000004D-5/
        DATA C /                                                          &
        .975083423708556D+0, -.240493938504146D-1,  .820452240880432D-3,  &
       -.434293081303427D-4,  .301844703403493D-5, -.025447331925082D-5,  &
        .002485835302051D-5, -.000273172013238D-5,  .000033084722797D-5,  &
       -.000004350549080D-5,  .000000614121457D-5, -.000000092236928D-5,  &
        .000000014635665D-5, -.000000002439278D-5,  .000000000424976D-5,  &
       -.000000000077084D-5,  .000000000014507D-5, -.000000000002824D-5,  &
        .000000000000567D-5, -.000000000000117D-5,  .000000000000025D-5,  &
       -.000000000000005D-5,  .000000000000001D-5/ 
      d=zero
      dd=zero
      if (abs(x) < three) then
!         CALCULATE VIA ERF
	 z=x/three
	 alpha=two-four*z*z
	 do j=na,0,-1
	    sv=d
	    d=-alpha*d-dd+a(J)
	    dd=sv
         end do
	 erfc=one-sqpi2*z*(d-dd)
      else
!      CALCULATE DIRECTLY
	 z=abs(three/x)
	 alpha=two-four*z*z
	 do j=nc,0,-1
	    sv=d
	    d=-alpha*d-dd+c(j)
	    dd=sv
         end do
	 if (x > zero) then 
	    erfc=half*exp(-x*x)/x*(d+half*alpha*dd)*sqpi2
         else
	    erfc=two-half*exp(-x*x)/(-x)*(d+half*alpha*dd)*sqpi2
         end if
      end if
    end function erfc
!
      subroutine inverse_matrix(A,n)
        implicit none
        integer(kind=i4_kind), intent(in)                  :: n
        real(kind=r8_kind), intent(inout),  dimension(n,n) :: a
        real(kind=r8_kind)                                 :: zero=0.0_r8_kind, &
                                                              one =1.0_r8_kind, &
                                                              d,tol,amax,t,swap,&
                                                              pivot
        integer(kind=i4_kind), allocatable, dimension(:,:) :: ipv
        integer(kind=i4_kind)                              :: alloc_stat,i,j,k, &
                                                              icolum,jcolum,    &
                                                              irow,jrow,l,l1,   &
                                                              nswap
        allocate(ipv(n,3), stat=alloc_stat)
        if(alloc_stat /= 0) &
             call error_handler("Inverse_matrix [1]: Allocation failed")
        D=1.D0
        TOL=1.0D-30
        DO J=1,N
           IPV(J,3)=0
        end DO

        DO I=1,N
           AMAX=0.D0
           DO J=1,N
              IF(IPV(J,3).EQ.1) cycle
              DO  K=1,N
                 IF(IPV(K,3).EQ.1) cycle
                 IF(AMAX.GE.ABS(A(J,K))) cycle
                 IROW=J
                 ICOLUM=K
                 AMAX=ABS(A(J,K))
              end DO
           end DO
           IF(AMAX.LE.TOL) then
              d=zero
              deallocate(ipv, stat=alloc_stat)
              if(alloc_stat /= 0) &
                   call error_handler("Inverse_matrix [2]: Deallocation failed")
              return
           end IF
           IPV(ICOLUM,3)=1
           IPV(I,1)=IROW
           IPV(I,2)=ICOLUM
           if(irow /= icolum) then
              DO  L=1,N
                 SWAP=A(IROW,L)
                 A(IROW,L)=A(ICOLUM,L)
                 A(ICOLUM,L)=SWAP
              end DO
           end if
           PIVOT=A(ICOLUM,ICOLUM)
           D=D*PIVOT
           A(ICOLUM,ICOLUM)=one
           DO L=1,N
              A(ICOLUM,L)=A(ICOLUM,L)/PIVOT
           end DO
           DO L1=1,N
              IF(L1.EQ.ICOLUM) cycle
              T=A(L1,ICOLUM)
              A(L1,ICOLUM)=0.D0
              DO L=1,N
                 A(L1,L)=A(L1,L)-A(ICOLUM,L)*T
              end DO
           end DO
        end DO

      NSWAP=0
      DO I=1,N
         L=N-I+1
         IF(IPV(L,1).EQ.IPV(L,2)) cycle
         JROW=IPV(L,1)
         JCOLUM=IPV(L,2)
         NSWAP=NSWAP+1
         DO K=1,N
            SWAP=A(K,JROW)
            A(K,JROW)=A(K,JCOLUM)
            A(K,JCOLUM)=SWAP
         end DO
      end DO
!c	if(mod(NSWAP,2).ne.0) D=-D
      D=D*((-1.d0)**NSWAP)
      deallocate(ipv, stat=alloc_stat)
      if(alloc_stat /= 0) &
           call error_handler("Inverse_matrix [2]: Deallocation failed")
    END subroutine inverse_matrix

    subroutine genmesh(rc,rsphere,ispher)
   
      implicit none
      real(kind=r8_kind), intent(in), dimension(3) ::  rc(3)         ! rc - centre of sphere 
      real(kind=r8_kind), intent(in)               ::  rsphere       ! radius 
      integer(kind=i4_kind), intent(in)            ::  ispher
      integer(kind=i4_kind)                        ::  l, m, alloc_stat
      real(kind=r8_kind)                           ::  a1(6,3),  a2(12,3)       
      real(kind=r8_kind)                           ::  a3(8,3),  b1(24,3)      
      real(kind=r8_kind)                           ::  b2(24,3), b3(24,3)       
      real(kind=r8_kind)                           ::  c1(24,3), k(3,3)        
      real(kind=r8_kind)                           ::  mk(3,3),  pk(3),qk(3)   
      real(kind=r8_kind)                           ::  lk(3,3)
      real(kind=r8_kind)                           ::  r1, r2                         
      real(kind=r8_kind)                           ::  asf(110,3)    ! weights fixed       
      integer(kind=4),dimension(3)                 ::  nspnt         ! 3 schemes are available
                                                              
      data lk /   0.00000000000000d0,    0.00000000000000d0,   0.00000000000000d0, &
                  0.30151134457800d0,    0.00000000000000d0,   0.00000000000000d0, &
                  0.18511563534500d0,    0.39568947305600d0,   0.69042104838200d0  /
      data mk /   0.00000000000000d0,    0.00000000000000d0,   0.00000000000000d0, &
                  0.90453403373300d0,    0.00000000000000d0,   0.00000000000000d0, &
                  0.96512403508700d0,    0.82876998125300d0,   0.21595729184600d0  /
      data pk /   0.88807383300000d0,    0.00000000000000d0,   0.87815891060400d0  /
      data qk /   0.45970084300000d0,    0.00000000000000d0,   0.47836902881200d0  /
      data asf/ 6*0.00952380952400d0,  8*0.03214285714300d0,24*0.02857142857100d0, &
               72*0.00000000000000d0,  6*0.01269841269800d0,12*0.02257495590800d0, &
                8*0.02109375000000d0, 24*0.02017333553800d0,60*0.00000000000000d0, &
                6*0.00382827049494d0,  8*0.00979373751249d0,24*0.00821173728319d0, &
               24*0.00959547133607d0, 24*0.00994281489118d0,24*0.00969499636166d0  /
      data nspnt / 38, 50, 110 /    
      allocate(screep_point(nspnt(ispher)), stat=alloc_stat)
      n_screep_points=0
      if(alloc_stat /= 0) &
           call error_handler("Genmesh [1]: Allocation failed")
      l=1
      call sfc(0.d0,0.d0,1.d0,0,0,1,a1,6,l)
      call sfc(0.d0,1.d0,0.d0,0,1,0,a1,6,l)
      call sfc(1.d0,0.d0,0.d0,1,0,0,a1,6,l)
      l=1
      r1=1.d0/sqrt(2.d0)
      call sfc(r1,r1,0.d0,1,1,0,a2,12,l)
      call sfc(r1,0.d0,r1,1,0,1,a2,12,l)
      call sfc(0.d0,r1,r1,0,1,1,a2,12,l)
      l=1
      r1=1.d0/sqrt(3.d0)
      call sfc(r1,r1,r1,1,1,1,a3,8,l)
      l=1
      r1=lk(1,ispher)
      r2=mk(1,ispher)
      call sfc(r1,r1,r2,1,1,1,b1,24,l)
      call sfc(r1,r2,r1,1,1,1,b1,24,l)
      call sfc(r2,r1,r1,1,1,1,b1,24,l)
      l=1
      r1=lk(2,ispher)
      r2=mk(2,ispher)
      call sfc(r1,r1,r2,1,1,1,b2,24,l)
      call sfc(r1,r2,r1,1,1,1,b2,24,l)
      call sfc(r2,r1,r1,1,1,1,b2,24,l)
      l=1
      r1=lk(3,ispher)
      r2=mk(3,ispher)
      call sfc(r1,r1,r2,1,1,1,b3,24,l)
      call sfc(r1,r2,r1,1,1,1,b3,24,l)
      call sfc(r2,r1,r1,1,1,1,b3,24,l)
      l=1
      r1=pk(ispher)
      r2=qk(ispher)
      call sfc(r1,r2,0.d0,1,1,0,c1,24,l)
      call sfc(r1,0.d0,r2,1,0,1,c1,24,l)
      call sfc(0.d0,r1,r2,0,1,1,c1,24,l)
      call sfc(r2,r1,0.d0,1,1,0,c1,24,l)
      call sfc(r2,0.d0,r1,1,0,1,c1,24,l)
      call sfc(0.d0,r2,r1,0,1,1,c1,24,l)

      do  m = 1, nspnt(ispher)
         screep_point(m)%square = 4*pi*asf(m,ispher)
      end do
write(output_unit,*) 'sum ', sum(screep_point(:)%square)
      call coor(rsphere,rc,a1, 6)     
      select case(ispher)
      case(1)
         call coor(rsphere,rc,a3, 8)
         call coor(rsphere,rc,c1,24)     
      case(2)
         call coor(rsphere,rc,a2,12)
         call coor(rsphere,rc,a3, 8)
         call coor(rsphere,rc,b1,24)
      case(3)
         call coor(rsphere,rc,a3, 8)
         call coor(rsphere,rc,b1,24)
         call coor(rsphere,rc,b2,24) 
         call coor(rsphere,rc,b3,24) 
         call coor(rsphere,rc,c1,24) 
      end select
                               
      if(n_screep_points .ne. nspnt(ispher) ) &
           call error_handler("Wrong number of points generated")
 
    end subroutine genmesh
  
    subroutine sfc(ax,ay,az,ix,iy,iz,a,ndim,l)
   
      implicit none
      integer(kind=i4_kind), intent(in)        :: ix,iy,iz,ndim
      integer(kind=i4_kind), intent(inout)     :: l
      integer(kind=i4_kind)                       isgnx,isgny,isgnz,nx,ny,nz
      integer(kind=i4_kind)                       ipx,ipy,ipz
      real(kind=r8_kind),    dimension(ndim,3) :: a
      real(kind=r8_kind),    intent(in)        :: ax,ay,az
      isgnx=1
      isgny=1
      isgnz=1
      nx=2
      ny=2
      nz=2
      if(ix.eq.0) nx=1
      if(iy.eq.0) ny=1
      if(iz.eq.0) nz=1
      do ipx=1,nx
         do ipy=1,ny
            do ipz=1,nz
               a(l,1)=ax*isgnx
               a(l,2)=ay*isgny
               a(l,3)=az*isgnz
               if(iz.ne.0) isgnz=-isgnz
               l=l+1
            end do
            if(iy.ne.0) isgny=-isgny
         end do
         if(ix.ne.0) isgnx=-isgnx
      end do
    
    end subroutine sfc

    subroutine coor(rad,rc,a,m)

      implicit none
      integer(kind=i4_kind)               :: i,m 
      real(kind=r8_kind), dimension(m,3)  :: a
      real(kind=r8_kind), dimension(3)    :: rc
      real(kind=r8_kind)                  :: rad

      do i=1,m
         n_screep_points=n_screep_points + 1
         screep_point(n_screep_points)%rc(1) = rc(1) + rad*a(i,1)
         screep_point(n_screep_points)%rc(2) = rc(2) + rad*a(i,2) 
         screep_point(n_screep_points)%rc(3) = rc(3) + rad*a(i,3)
      end do
    end subroutine coor
!
!   Recursive calculations for Normalized Real Solid Harmonics / r**l
!
      subroutine sh(x,y,z,lmax,c)

        implicit none
        real(kind=r8_kind), dimension(:,:) :: c ! ylm(1:lmax+1, 1:2*lmax+1)
        real(kind=r8_kind), intent(in)     :: x,y,z
        real(kind=r8_kind)                 :: r,sqr
        integer(kind=i4_kind),intent(in)   :: lmax
        integer(i4_kind)                   :: l,m,im,i

        c(1,1) = 1.d0 
        if(lmax==0) return
        r = sqrt(x**2 +y**2 + z**2)   
        if(r==0.d0) then
           c(2:lmax+1,:)=0.d0
           return
        end if
        c(2,1) = z / r
        c(2,2) = x / r
        c(2,3) = y / r
        if(lmax==1) return 
        do l=3,lmax+1
           m=1
           c(l,m)=p(l-1,z/r)
           do im=1,l-3
              do i=1,2
                 m=m+1   
                 c(l,m) = 1.d0/sqrt( float( (l-1)**2-(m/2)**2) ) *      &
                      ( (2*l-3) * z/r * c(l-1,m) -                      &
                      sqrt( float( (l-2)**2-(m/2)**2 ) ) * c(l-2,m) )  
              end do
           end do  ! im
           im=l-2
           do i=1,2
              m=m+1
              c(l,m) = 1.d0/sqrt( float( (l-1)**2-(m/2)**2) ) *          & 
                   (2*l-3) * z/r * c(l-1,m)         
           end do
           im=l-1
           sqr=sqrt( float(2*l-3)/float(2*l-2) )
           m=m+1
           c(l,m) = sqr * ( x*c(l-1,2*l-4) - y*c(l-1,2*l-3) ) / r                  
           m=m+1
           c(l,m) = sqr * ( y*c(l-1,2*l-4) + x*c(l-1,2*l-3) ) / r       
        end do  ! l
        return
      contains
  
        function p(n,t)
          implicit none
          real(kind=r8_kind)    :: p,t,pk,pk1,pk2
          integer(kind=i4_kind) :: n,k

          select case(n)
          case(0)
             p=1.d0
             return
          case(1)
             p=t
             return
          end select
          pk2=1.d0
          pk1=t
          do  k=2,n
             pk=( (2*k-1)*t*pk1 - (k-1)*pk2 ) / k
             pk2=pk1
             pk1=pk
          end do
          p=pk
          return
        end function p

      end subroutine sh

      subroutine project(vs,lmax,rc,c)      
   
        implicit none
        real(kind=r8_kind), dimension(:,:)           :: c 
        real(kind=r8_kind), intent(in)               :: vs(:)   
        real(kind=r8_kind), intent(in), dimension(3) :: rc
        integer(kind=i4_kind), intent(in)            :: lmax
        integer(kind=i4_kind)                        :: l,m,nlm,ip
 
        vlm(:)=0.d0
        do ip=1,n_screep_points
           call sh(screep_point(ip)%rc(1)-rc(1), &
                   screep_point(ip)%rc(2)-rc(2), &
                   screep_point(ip)%rc(3)-rc(3),lmax,c)
           nlm=0
           do l=1,lmax+1
              do m=1,2*l-1
                 nlm=nlm + 1
                 vlm(nlm) = vlm(nlm) + vs(ip) * c(l,m) * screep_point(ip)%square
              end do
           end do
        end do

        write(output_unit,'(/30(1h*), " Y_lm amlitudes ",20(1h*)/&
                          &t10,"N", t15,"l", t20,"m", t35,"Vlm")' )
        nlm=0
        do l=1,lmax+1
           do m=1,2*l-1
              nlm=nlm+1
              write(output_unit,'( t8,i3,t14,i2,t19,i2,t27,f19.15,3x)') & 
                                nlm,l,m,vlm(nlm)
           end do
        end do
        write(output_unit,'(/67(1h*)/)')
      
      end subroutine project
!
!     Restore potential in (x,y,z)-point inside a sphere
!
      subroutine rest(x,y,z,rs,rc,lmax,c,vrest)

        implicit none 
        real(kind=r8_kind), dimension(:,:)            :: c 
        real(kind=r8_kind), intent(in), dimension(3)  :: rc
        real(kind=r8_kind), intent(in)                :: x,y,z,rs
        real(kind=r8_kind)                            :: r,rrl
        real(kind=r8_kind), intent(out)               :: vrest
        integer(kind=i4_kind), intent(in)             :: lmax
        integer(kind=i4_kind)                         :: lm,l,m

        call sh(x-rc(1),y-rc(2),z-rc(3),lmax,c)    
        vrest = 0.d0
        rrl=1.d0
        r=sqrt( (x-rc(1))**2 + (y-rc(2))**2 + (z-rc(3))**2 )
        lm=0
        do l=1,lmax+1
           if(l .gt. 1) rrl = rrl * (r/rs) 
           do m=1,2*l-1
              lm=lm+1
              vrest = vrest + rrl * float(2*l-1)/(4*pi) * c(l,m)*vlm(lm)
           end do
        end do

      end subroutine rest

    end module ewald_module
!
!
!
      program main
        use ewald_module

        call ewald

      end program main







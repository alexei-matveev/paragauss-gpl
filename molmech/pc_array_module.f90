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
module pc_array_module

  !------------ Modules used -----------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use species_module

  implicit none
  private
  save
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of public constants and variables -----
  !------------ public functions and subroutines ---------------------
  public read_pc_array_options, read_embed_cluster, gen_pc_array
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of private constants and variables ----
  type screep_point_type
     real(r8_kind), dimension(3) :: rc
     real(r8_kind)               :: V_dif,q,square
  end type screep_point_type

  type atom_type
     real(r8_kind), dimension(3) :: s    ! shell  coordinate
     real(r8_kind), dimension(3) :: c    ! core  coordinate
     real(r8_kind),  dimension(3) :: rc   ! cartesian coordinate
     real(r8_kind), dimension(3) :: rl   ! lattice   coordinate
     real(r8_kind)               :: q    ! charge
     integer(i4_kind)            :: qt   ! charge type
     real(r8_kind)               :: an   ! real for charge type
     character(len=len_name)     :: name ! atom symbol
  end type atom_type

  type cluster_type
     real(r8_kind), dimension(3) :: rc   ! cartesian coordinate
     real(r8_kind), dimension(3) :: rl   ! lattice   coordinate
     real(r8_kind)               :: V_pc ! PC-array potential on ato
     real(r8_kind)               :: V_ew ! madelung potential on atom
     character(len=2)                 :: name ! atom symbol
  end type cluster_type

  type lattice_type
     real(r8_kind),    pointer, dimension(:,:) :: v_cartesian
     real(r8_kind),    pointer, dimension(:)   :: v_mod
     integer(i4_kind), pointer, dimension(:,:) :: v_vector
     integer(i4_kind)                          :: num_vector
     integer(i4_kind), pointer, dimension(:)   :: num_in_sphere
     integer(i4_kind)                          :: n_spheres
  end type lattice_type

  type point_charge
     real(r8_kind), dimension(3) :: rc
     real(r8_kind)               :: q
     character(len=2)                 :: name
  end type point_charge

  type(screep_point_type), pointer, dimension(:)     :: screep_point
  type(atom_type),         pointer, dimension(:)     :: atomew
  type(cluster_type),      pointer, dimension(:)     :: cluster
  type(point_charge),      pointer, dimension(:)     :: pc
  real(r8_kind),  allocatable, dimension(:,:)   :: vertex,rl
  real(r8_kind),  allocatable, dimension(:)     :: vlm
  real(r8_kind)                  :: V_cell, V_BZ, eta
  real(r8_kind), dimension(3)    :: Rc_screepAU,Rc_screepANGS
  real(r8_kind), dimension(3,3)  :: a, b
  integer(i4_kind)               :: n_triangle, n_atoms_in_cell, n_cluster_atoms
  integer(i4_kind), parameter    :: input_unit=5, output_unit=6 ,gxc_unit=7
  integer(i4_kind)               :: N_ewald_pc, n_div
  integer(i4_kind)               :: i_scheme,lmax
  character(len=100)                  :: work_dir, input_file, output_file
  character(len=100)                  :: screep_procedure
  logical                             :: axes_c3, axes_c4, input_file_exist
  logical                             :: symmetry_c3v=.false.

  type(lattice_type)                  :: lattice
  real(r8_kind)                       :: R_max_fragment
  logical, parameter                  :: print_vectors=.true.
  logical, parameter                  :: not_print_vectors=.false.
  logical, parameter                  :: full_print=.true.
  logical                             :: make_crystal_fragment, make_ewald_array ! tasks
!  logical:: gxc_file_exist
!  logical:: gxcv_file_exist

  logical, dimension(:,:), allocatable:: type_number

  real(r8_kind) :: r_max_pc
  real(r8_kind) :: r_screep
  integer(i4_kind) :: n_screep_points
  real(r8_kind) :: rc_screep(3)
  character(len=3) :: axes
  real(r8_kind) :: direction(3)

  real(r8_kind) :: df_r_max_pc=17.5_r8_kind
  real(r8_kind) :: df_r_screep=five
  integer(i4_kind) :: df_n_screep_points=128
  real(r8_kind) :: df_rc_screep(3)=zero
  character(len=3) :: df_axes='C4'
  real(r8_kind) :: df_direction(3)=(/zero,zero,one/)

  namelist/pc_array_data/ r_max_pc,r_screep,rc_screep,n_screep_points,axes,direction


  !------------ Subroutines -----------------------------------------
contains
  !******************************************************************
  function read_pc_array_options()

    logical :: read_pc_array_options

    integer(i4_kind) :: i

    axes_c3 =.false.
    axes_c4 =.false.
    eta=0.17_r8_kind

    r_max_pc=df_r_max_pc
    r_screep=df_r_screep
    n_screep_points=df_n_screep_points
    rc_screep=df_rc_screep
    axes=df_axes

    call  go_to_first_input_line
    read_pc_array_options=find_namelist("&PC_ARRAY_DATA",i)
    if(read_pc_array_options) read(input_device,nml=pc_array_data, end=100, err=200)
100 call upcase(axes)
    select case(trim(axes))
    case("C3")
       axes_c3=.true.
    case("C3V")
       axes_c3=.true.
       symmetry_c3v=.true.
    case("C4")
       axes_c4=.true.
    end select
    n_triangle=8

    n_div = int (sqrt (float (n_screep_points / n_triangle))) - 1
    n_screep_points=n_triangle*(n_div+1)**2
    direction(:)=direction(:)/ &
         sqrt(dot_product(direction,direction))

    return

200 call input_nm_error(0,"PC_ARRAY_DATA")

  end function  read_pc_array_options
  !******************************************************************

  !******************************************************************
  function read_embed_cluster()

    logical :: read_embed_cluster

    logical :: start_box,end_box
    integer(i4_kind) :: i0,i1,i_atom,alloc_stat
    real(r8_kind) :: tmp_r(3)
    character(len=6) :: number

    i0=0; i1=0

    start_box=find_word("&EMBED_CLUSTER",i0)
    end_box=find_word("/EMBED_CLUSTER",i1)

    if(.not. start_box .or. .not. end_box) &
         call error_handler("MolMech: Not terminated EMBED_CLUSTER box of data")

    n_cluster_atoms=i1-i0-1
    if(n_cluster_atoms <= 0) n_cluster_atoms=0

    if(.not. start_box .and. .not. end_box) then
       read_embed_cluster=.false.
       return
    end if

    allocate(cluster(n_cluster_atoms), stat=alloc_stat)
    if(alloc_stat /= 0) &
         call error_handler("MolMech: failed CLUSTER allocation")
    start_box=find_word("&EMBED_CLUSTER",i0)
    do i_atom=1,n_cluster_atoms
       read(input_device,*, err=300)  cluster(i_atom)%name,tmp_r
       cluster(i_atom)%rc=(tmp_r-Rc_screep)
    end do
    do i_atom=1,n_cluster_atoms
       if(sqrt(dot_product((cluster(i_atom)%rc), &
            (cluster(i_atom)%rc)) )>=r_screep) then
          call error_handler("MolMech: Cluster atom is out of sphere")
       end if
    end do

    read_embed_cluster=.true.
    return

300 write(number,'(i4)') i_atom
    call error_handler("MolMech: Coordinates of "//trim(number)// &
         " cluster atom has been read in with error")

  end function read_embed_cluster
  !******************************************************************

  !******************************************************************
  subroutine init_data()

    integer(i4_kind) :: i,j,k,alloc_stat,i_atom,typew
    real(r8_kind) :: tmp_r(3)

    eta=0.17_r8_kind

    a(:,1)=vect%v1
    a(:,2)=vect%v2
    a(:,3)=vect%v3

    V_cell=dot_product(a(:,1),cross_product(a(:,2),a(:,3)))
    !Reciprocal lattice vectors
    do i=1,3
       j=i+1
       if(j.gt.3) j=j-3
       k=j+1
       if(k.gt.3) k=k-3
       b(:,i)=2*pi*cross_product(a(:,j),a(:,k) )/V_cell
    end do
    V_BZ=8*pi**3/V_cell

    do i_atom=1,n_cluster_atoms
       cluster(i_atom)%rl=matmul(cluster(i_atom)%rc,b)/(two*pi)
    end do

    n_atoms_in_cell=n_species

    allocate(atomew(n_atoms_in_cell),stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: failed ATONEW allocation")

    do i_atom=1,n_atoms_in_cell
       tmp_r=atoms_cart(i_atom)%r
       typew=atoms_cart(i_atom)%type
       atomew(i_atom)%name=atoms(typew)%name
       atomew(i_atom)%q=atoms(typew)%charge
       atomew(i_atom)%rc=(tmp_r-Rc_screep)
       atomew(i_atom)%rl= &
            matmul(tmp_r-Rc_screep,b)/(two*pi)
    end do
    Rc_screepAU=Rc_screep

  end subroutine init_data
  !******************************************************************

  !******************************************************************
  subroutine gen_pc_array()

!!$    character(len=100)                  :: gxc_file
!!$    character(len=100)                  :: gxcv_file
    real(r8_kind), pointer,dimension(:)        :: V_ew
    real(r8_kind), allocatable, dimension(:,:) :: A
    real(r8_kind)                              :: R
    real(r8_kind)                              :: cl_V
    integer(i4_kind)                           :: i_atom,i_point,alloc_stat,i,j

!!$    gxc_file=trim(work_dir)//"/"//"gxcfile"
!!$    gxcv_file=trim(work_dir)//"/"//"gxcvfile"

!!$    inquire(file=gxc_file, exist=  gxc_file_exist)
!!$    inquire(file=gxcv_file, exist=  gxcv_file_exist)
!!$    if(gxcv_file_exist)  then
!!$       open(gxc_unit,file=trim(gxcv_file),status="old")
!!$    else if(gxc_file_exist) then
!!$       open(gxc_unit,file=trim(gxc_file),status="old")
!!$    endif

    call init_data()

    call pc_generate()

    allocate(rl(n_cluster_atoms,3), stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed RL allocation(1)")
    do i_atom=1,n_cluster_atoms
       rl(i_atom,:)=cluster(i_atom)%rl(:)
    end do
    V_ew => cluster(:)%V_ew
print*,'Madelung'
    call Madelung(rl,V_ew,n_cluster_atoms)

    deallocate(rl, stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed RL deallocation(1)")


print*,'screep_generate'
    call screep_generate()

    allocate(rl(n_screep_points,3), stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed RL allocation(2)")
    do i_point=1,n_screep_points
       rl(i_point,:)=matmul(screep_point(i_point)%rc,b)/(two*pi)
    end do
    V_ew => screep_point(:)%V_dif
print*,'Madelung'
    call Madelung(rl,V_ew,n_screep_points)
    deallocate(rl, stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed RL deallocation(2)")

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

    allocate(A(n_screep_points, n_screep_points), stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed A deallocation")

    do i=1,n_screep_points
       A(i,i)=1.07_r8_kind*sqrt(four*acos(-one)/screep_point(i)%square)
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

    deallocate(A,stat=alloc_stat)
    if(alloc_stat /= 0) call error_handler("MolMech: PC_ARRAY_MODULE: failed A deallocation")

!!$            write(output_unit,"(/10X, 'Total charge on sphere :', f9.5)") &
!!$                    sum(screep_point(:)%q)
!!$            do i=1,n_screep_points
!!$               write(output_unit,*) i, 'q=', screep_point(i)%q
!!$            end do

    do i_atom=1,n_cluster_atoms
       cl_V=cluster(i_atom)%V_pc
       do j=1,n_screep_points
          R=sqrt(dot_product(screep_point(j)%rc-cluster(i_atom)%rc, &
               screep_point(j)%rc-cluster(i_atom)%rc  ) )
          if(R <= 1.0e-10_r8_kind) R=zero
          if(R /= zero) cl_V=cl_V + screep_point(j)%q/R
       end do
!!$       format="(10x,i3,3x,3(f15.7,3x))"
!!$               write(output_unit,format) i_atom,cluster(i_atom)%V_ew*27.2113961d0, cl_V*27.2113961d0, &
!!$                                                abs(cluster(i_atom)%V_ew-cl_V)*27.2113961d0
       write(*,*) i_atom,cluster(i_atom)%V_ew*14.3997584_r8_kind, cl_V*14.3997584_r8_kind, &
            abs(cluster(i_atom)%V_ew-cl_V)*14.3997584_r8_kind
    end do
!++++++++++++++++++++++++++++ Store embedded cluster ++++++++++++++++++++++++++
!!$            open(20, file="ew_pc.xyz", status="unknown")
!!$            open(30, file="ewald.pcr", status="unknown")
!!$            write(20,*) n_screep_points+n_ewald_pc
!!$            write(30,*) n_screep_points+n_ewald_pc
!!$            write(20,*) "comment line"
!!$            i_at: do i_atom=1,n_ewald_pc
!!$               do i=1,n_cluster_atoms
!!$                  if( all(abs(pc(i_atom)%rc-cluster(i)%rc)<=1.d-10) ) then
!!$                     write(20,"(' Si',4x,3(f12.7,3x),15x,3x,i4)") &
!!$                             0.529177_r8_kind* (cluster(i)%rc(:)+Rc_screepAU), i
!!$                    write(30,"(4f15.8,i3,4i2,6x,f4.2)") &
!!$                    pc(i_atom)%rc(:)+Rc_screepAU, pc(i_atom)%q,1,0,0,1,0, 0.0
!!$                     cycle i_at
!!$                  end if
!!$               end do
!!$               write(20,"(a4,4x,3(f12.7,3x),5x,f10.7,3x,i4)") &
!!$                    pc(i_atom)%name, &
!!$                    pc(i_atom)%rc(:)+Rc_screepAU, pc(i_atom)%q, i_atom
!!$                write(30,"(4f15.8,i3,4i2,6x,f4.2)") &
!!$                    pc(i_atom)%rc(:)+Rc_screepAU, pc(i_atom)%q,1,0,0,1,0, 0.0
!!$            end do i_at
!!$
!!$            do i=1,n_screep_points
!!$               write(20,"('H',3(5x,f12.7))" ) &
!!$                    RC_screepAU+screep_point(i)%rc(:)
!!$              write(30,"(4f15.8,i3,4i2,6x,f4.2)") &
!!$              screep_point(i)%rc(:)+Rc_screepAU, &
!!$              screep_point(i)%q,1,0,0,1,0, 0.0
!!$            end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    call output_pc_array_data()

    nullify(V_ew)

    call shudown_all_arrays()

  end subroutine gen_pc_array
  !******************************************************************

  !******************************************************************
  subroutine output_pc_array_data()

    integer(i4_kind) :: i,j
    real(r8_kind) :: cl_V,R

    write(output_device,'(80("*"))')
    write(output_device,'(25x,a29,26x)') 'RESULT OF PC ARRAY GENERATION'
    write(output_device,'(/)')
    write(output_device,'(a17)') 'Embedded cluster:'
    do i=1,n_cluster_atoms
       write(output_device,'(i4,a6,3f15.7)') i,cluster(i)%name,cluster(i)%rc+Rc_screep
    end do
    write(output_device,'(/)')
    write(output_device,'(a14)') 'SCREEP sphere:'
    write(output_device,'(a8,3f15.7)') 'Center: ',Rc_screep
    write(output_device,'(a8,f15.7)') 'Radius: ',R_screep
    write(output_device,'(a18,i6)') 'Number of points: ',n_screep_points
    write(output_device,'(a10,a4)') 'Symmetry: ',axes
    write(output_device,'(/)')
    write(output_device,'(a9)') 'PC array:'
    write(output_device,'(a25,i7)') 'Number of point charges: ',n_screep_points+n_ewald_pc
    write(output_device,'(/)')
    write(output_device,'(a22)') 'Result of fitting (V):'
    write(output_device,'(a12,2x,a18,2x,a16,2x,2x,a10)') &
         'Cluster atom','Madelung potential','Fitted potential','Difference'
    write(output_device,'(a64)') &
         '-----------------------------------------------------------------'
    do i=1,n_cluster_atoms
       cl_V=cluster(i)%V_pc
       do j=1,n_screep_points
          R=sqrt(dot_product(screep_point(j)%rc-cluster(i)%rc, &
               screep_point(j)%rc-cluster(i)%rc))
          if(R <= 1.0e-10_r8_kind) R=zero
          if(R /= zero) cl_V=cl_V + screep_point(j)%q/R
       end do
       write(output_device,'(4x,i4,4x,1x,f15.7,2x,1x,f15.7,3x,f15.7)') &
            i,cluster(i)%V_ew*14.3997584_r8_kind, &
            cl_V*14.3997584_r8_kind, &
            abs(cluster(i)%V_ew-cl_V)*14.3997584_r8_kind
    end do



    write(output_device,'(80("*"))')

  end subroutine output_pc_array_data
  !******************************************************************

  !******************************************************************
  subroutine shudown_all_arrays

    integer(i4_kind) :: status

    deallocate(cluster,atomew,screep_point,pc, stat=status)
    if(status /= 0) call error_handler("MolMech: PC_ARRAY_MODULE: failed CLUSTER deallocation")

  end subroutine shudown_all_arrays
  !******************************************************************

  !******************************************************************
  subroutine screep_generate

    integer(i4_kind), allocatable, dimension(:,:) :: index
    integer(i4_kind)                              :: alloc_stat,i_point
    integer(i4_kind)                              :: N_vertex,i_triangle
    integer(i4_kind)                              :: i_div, step_index
    integer(i4_kind)                              :: n_step, i_vertex
    real(r8_kind), dimension(3)                   :: v_1, v_2
    real(r8_kind), dimension(3)                   :: p_0,p_1,p_2,p_3
    real(r8_kind)                                 :: tmp_axes(3), tmp_vector(3)
    real(r8_kind), dimension(3)                   :: rot_axes
    real(r8_kind), dimension(3,3)                 :: u
    real(r8_kind), allocatable, dimension(:,:)    :: mv1, mv2, mv3
    real(r8_kind)                                 :: scalar,delta,square,s_sum
    real(r8_kind)                                 :: sp_square, sp_sum
    real(r8_kind)                                 :: mean, mean2
!!$    character(len=100)                                 :: format,screep_file

!!$         screep_file=trim(work_dir)//"/"//"screep.xyz"
!!$         screep_unit=8      ! to be changed by "open_get ..."
!!$         open(screep_unit,file=trim(screep_file), status="unknown")

    allocate(screep_point(n_screep_points), stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed SCREEP_POINT allocation")

    allocate(index(n_triangle,3),stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed INDEX allocation")
    index=reshape((/1,2,3,4,1,2,3,4,2,3,4,1, &
         2,3,4,1,5,5,5,5,6,6,6,6/),(/n_triangle,3/))

    N_vertex=(maxval(index))
    allocate(vertex(3,N_vertex), stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed VERTEX allocation")
    !         C_4 || z-axes (default, will be changed depending on axes chosen
    vertex=reshape((/ 1.0_r8_kind, 0.0_r8_kind, 0.0_r8_kind,    &
         0.0_r8_kind, 1.0_r8_kind, 0.0_r8_kind,    &
         -1.0_r8_kind, 0.0_r8_kind, 0.0_r8_kind,    &
         0.0_r8_kind,-1.0_r8_kind, 0.0_r8_kind,    &
         0.0_r8_kind, 0.0_r8_kind, 1.0_r8_kind,    &
         0.0_r8_kind, 0.0_r8_kind,-1.0_r8_kind /), &
         (/3,N_vertex/) )

    if(.not.axes_c3 .and. .not.axes_c4 ) &
         call error_handler("MolMech: PC_ARRAY_MODULE: axes is not defined")

    if(axes_c4) tmp_axes=(/0.0_r8_kind, 0.0_r8_kind, 1.0_r8_kind/)
    if(axes_c3) then
       tmp_axes(:)=(vertex(:,1)+vertex(:,2)+vertex(:,5))/3.0_r8_kind
       scalar=sqrt( dot_product(tmp_axes,tmp_axes) )
       tmp_axes(:)=tmp_axes(:)/scalar
    end if

    delta=acos( dot_product(direction,tmp_axes)/ &
         (sqrt(dot_product(direction,direction))* &
         sqrt(dot_product(tmp_axes,tmp_axes))))

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

    !
    !Fill the triangle-faces
    !
    allocate(mv1(3,n_screep_points),mv2(3,n_screep_points),mv3(3,n_screep_points), &
         stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed MV1, MV2, MV3 allocation")

    i_point=0
    triangle:  do i_triangle=1,n_triangle
       v_1(:)=(vertex(:,index(i_triangle,2))-vertex(:,index(i_triangle,1)))/&
            (n_div+1)
       v_2(:)=(vertex(:,index(i_triangle,3))-vertex(:,index(i_triangle,1)))/&
            (n_div+1)

       do i_div=1,n_div+1
          p_0(:)=vertex(:,index(i_triangle,1))+add_vector(0,v_1,i_div-1,v_2,3)
          step_index = -1
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
       s_sum=s_sum+square
       sp_sum = sp_sum + sp_square
       mean=mean + sp_square
       mean2=mean2 + sp_square**2
!!$           format="(10X/i4,5x,'S=',f10.5,5x,'Sph=',f10.5)"
!!$
!!$           write(output_unit,format) i_point,square,sp_square
!-----------------
    end do

!!$        write(output_unit,*) s_sum, sp_sum, four*pi*R_screep**2
!!$        write(output_unit,*) "Dispersion : ", &
!!$                           sqrt(mean2/n_screep_points - (mean/n_screep_points)**2)

    deallocate(mv1,mv2,mv3,vertex,index, stat=alloc_stat)
    if(alloc_stat.ne.0) &
         call error_handler("MolMech: PC_ARRAY_MODULE: failed MV1, MV2, MV3 deallocation")

!!$!                    Save generated screep_points in MOLDEN ".xyz" format
!!$        write(screep_unit,*) n_screep_points
!!$        write(screep_unit,*) "comment line"
!!$        do i_point=1,n_screep_points
!!$           format="('H',3(5x,f12.7))"
!!$         write(screep_unit,format) &
!!$              screep_point(i_point)%rc(:)+Rc_screepAU
!!$        end do
!!$        close(screep_unit, status='keep')

  end subroutine screep_generate
  !******************************************************************

  !******************************************************************
  function add_vector(n_1,vec_1,n_2,vec_2,dim)

    integer(i4_kind), intent(in)            :: n_1,n_2,dim
    real(r8_kind), dimension(3), intent(in) :: vec_1,vec_2
    real(r8_kind), dimension(dim)             :: add_vector

    add_vector(:)=n_1*vec_1(:) + n_2*vec_2(:)

  end function add_vector
  !******************************************************************

  !******************************************************************
  function cross_product(v_1,v_2)

    real(r8_kind)               :: cross_product(3)
    real(r8_kind), dimension(3) :: v_1, v_2

    cross_product(1)=v_1(2)*v_2(3)-v_1(3)*v_2(2)
    cross_product(2)=v_1(3)*v_2(1)-v_1(1)*v_2(3)
    cross_product(3)=v_1(1)*v_2(2)-v_1(2)*v_2(1)

  end function cross_product
  !******************************************************************

  !******************************************************************
  subroutine rot(u,rot_axes,d)
    !======================================================================!
    !                  rotation around axes (nx,ny,nz) on d                !
    !======================================================================!
    real(r8_kind),intent(inout), dimension(3)   :: rot_axes(3)
    real(r8_kind),intent(inout), dimension(3,3) :: u
    real(r8_kind),intent(in)                    :: d
    real(r8_kind)                               :: cosd,sind,nx,ny,nz

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
  !******************************************************************

  !******************************************************************
  function sph_square(v_1,v_2,v_3,R)

    real(r8_kind)                           :: sph_square
    real(r8_kind), intent(in), dimension(3) :: v_1,v_2,v_3
    real(r8_kind), intent(in)               :: R
    real(r8_kind)                           :: alpha, beta, gamma

!!$    pi=acos(-1.0_r8_kind)
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
  !******************************************************************

  !******************************************************************
  subroutine generate_lattice(space,mu_1,mu_2,mu_3,r_max,print_option,full_print)
    !===================================================================================!
    !                                                                                   !
    !       Purpose: subroutine generates the set of lattice vectors,                   !
    !                both in real and reciprocal space;                                 !
    !                In case of rec.space, the variables mu1, mu2, mu3                  !
    !                define lattice coordinates of k-vector (or k-point)                !
    !                If option "full_print" is included the full output will            !
    !                be generated (i.e. vector, cartesian coordiates and                !
    !                modules of lattice vetors)                                         !
    !===================================================================================!
    real(r8_kind),      intent(in)      :: mu_1, mu_2, mu_3
    real(r8_kind),      intent(in)      :: r_max
    logical,           intent(in)      :: print_option
    logical, optional, intent(in)      :: full_print
    character(len=*),  intent(in)      :: space ! 'direct' or 'reciprocal'

    !                              local variables

    real(r8_kind)                                 :: aux_vector(3),R(3),h(6)
    real (r8_kind), dimension(3) :: v_temp_cartesian
    integer (i4_kind), dimension(3) :: v_temp_vector
    real(r8_kind), dimension(3,3)                 :: a_lat,v
    real(r8_kind)                                 :: a_min, mod_temp
    integer(i4_kind)                              :: i,k,l
    integer(i4_kind)                              :: num_slab, num_in_slab
    integer(i4_kind)                              :: n_1, n_2, n_3
    integer(i4_kind)                              :: i_cycle, num, j
    integer(i4_kind)                              :: i_vector, j_vector
    logical                                            :: empty_cycle
!!$    character(len=32) print_format

    select case(space)
    case('direct')
       a_lat = a
    case('reciprocal')
       a_lat = b
    end select

    !      First, empty cycle to obtain the number of vectors (num_vector)
    !                          Inscribe a "R_max sphere"
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

!!$    print_format="(10x, 'Number of vectors : ', i5/10x, 'Number of shells :',i5 )"
!!$    write(output_unit,*) lattice%num_vector, lattice%n_spheres
!!$    do i_print=1,(lattice%n_spheres-1)/10+1
!!$       print_format="(10x,10i5)"
!!$       write(output_unit,print_format) &
!!$            (lattice%num_in_sphere(m),m=((i_print-1)*10+1),&
!!$            min(i_print*10,lattice%n_spheres))
!!$    end do
!!$
!!$    if(print_option) then
!!$       if( present(full_print) ) then
!!$          write(6,'(/20x,"Coordinates in vector units"/)')
!!$          do i_print=1,lattice%num_vector
!!$             i_first=(i_print-1)*5 + 1
!!$             i_second=i_print*5
!!$             if(i_second > lattice%num_vector) i_second=lattice%num_vector
!!$             print_format='(3x, 5("(",3i4,")",1x), i3)'
!!$             write(6,print_format) ( (lattice%v_vector(m,j),j=1,3), &
!!$                  m=i_first,i_second),i_second
!!$             if(i_second==lattice%num_vector) exit
!!$          end do
!!$
!!$          write(6,'(/20x,"Cartesian coordinates"/)')
!!$          do i_print=1,lattice%num_vector
!!$             i_first=(i_print-1)*5 + 1
!!$             i_second=i_print*5
!!$             if(i_second > lattice%num_vector) i_second=lattice%num_vector
!!$             print_format='( 3x, 5("(",3f9.4,")",2x), / )'
!!$             write(6,print_format) ( (lattice%v_cartesian(m,j),j=1,3), &
!!$                  m=i_first,i_second)
!!$             if(i_second==lattice%num_vector) exit
!!$          end do
!!$       end if
!!$
!!$       write(6,'(/20x, "Vector modules"/)')
!!$       do i_print=1,lattice%num_vector
!!$          i_first=(i_print-1)*5 + 1
!!$          i_second=i_print*5
!!$          if(i_second > lattice%num_vector) i_second=lattice%num_vector
!!$          print_format='(3x, 5g12.5, i5)'
!!$          write(6,print_format) (lattice%v_mod(m), m=i_first,i_second), i_second
!!$          if(i_second==lattice%num_vector) exit
!!$       end do
!!$    end if ! print_option

   end subroutine generate_lattice
  !******************************************************************

  !******************************************************************
   subroutine pc_generate

     type(atom_type),    pointer, dimension(:)     :: atom_wz,tmpe
     type(point_charge), pointer, dimension(:)     :: tmp_pc
     real(r8_kind)                            :: mu_1,mu_2,mu_3,q_cell,q_tot,r_max,R
     real(r8_kind),dimension(3)               :: p_cell,rs
     integer(i4_kind)                         :: N_wz,i_atom,n_equiv,i,j,k
     integer(i4_kind)                         :: alloc_stat

!!$        scell_file=trim(work_dir)//"/"//"cell.xyz"
!!$        scell_unit=7                                    ! to be changed by "openget ..."
!!$        open(scell_unit, file=trim(scell_file), status="unknown")
!!$        pc_xyz_file=trim(work_dir)//"/"//"Ewald_PC.xyz"
!!$        pc_xyz_unit=10                                  ! to be changed by "openget ..."
!!$        open(pc_xyz_unit, file=trim(pc_xyz_file), status="unknown")
!                                  Form symmetric cell
     N_wz=0
     do i_atom=1,n_atoms_in_cell
        mu_1=atomew(i_atom)%rl(1)
        mu_2=atomew(i_atom)%rl(2)
        mu_3=atomew(i_atom)%rl(3)
        r_max=sqrt( dot_product(atomew(i_atom)%rc,atomew(i_atom)%rc) )
        if(r_max==0.0_r8_kind) r_max=1.0_r8_kind
        do
           call generate_lattice("direct",mu_1,mu_2,mu_3,r_max,print_vectors)
           if(lattice%n_spheres >= 2) exit
           r_max=2*r_max
        end do        !          First shell only
!!$        write(output_unit,*) 'atom ',i_atom, 'n_equiv= ',lattice%num_in_sphere(1)
        n_equiv=lattice%num_in_sphere(1)
        if(i_atom > 1) then
           allocate(tmpe(N_wz), stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("MolMech: PC_ARRAY_MODULE: failed TMPE allocation")
           tmpe=atom_wz
           deallocate(atom_wz, stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("MolMech: PC_ARRAY_MODULE: failed ATOM_WZ deallocation")
           allocate(atom_wz(N_wz+n_equiv), stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("MolMech: PC_ARRAY_MODULE: failed ATOM_WZ allocation")
           do i=1,N_wz
              atom_wz(i)=tmpe(i)
           end do
           deallocate(tmpe, stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("MolMech: PC_ARRAY_MODULE: failed TMPE deallocation")
        else
           allocate(atom_wz(n_equiv), stat=alloc_stat)
           if(alloc_stat.ne.0) &
                call error_handler("MolMech: PC_ARRAY_MODULE: failed ATOM_WZ allocation(1)")
        end if
        do i=1,n_equiv
           atom_wz(N_wz+i)%rc(:)=lattice%v_cartesian(i,:)
           atom_wz(N_wz+i)%rl(:)=atomew(i_atom)%rl(:)+lattice%v_vector(i,:)
           atom_wz(N_wz+i)%q=atomew(i_atom)%q / n_equiv
           atom_wz(N_wz+i)%name=atomew(i_atom)%name
        end do
        N_wz=N_wz+n_equiv
     end do ! i_atom
!                                Control cell-charge and store sym.cell
     q_cell=0.0_r8_kind
     p_cell=0.0_r8_kind
!!$        write(output_unit,"(/37(1h#),' Symmetric cell ',37(1h#)  )" )
!!$        write(output_unit,"(10X,'Symmetric cell is stored in : ',a)") scell_file
!!$        write(scell_unit,"(i7/'comment line')" ) N_wz
!!$        write(output_unit,"(/30X,'X',17x,'Y',20x,'Z',17x,'Q' )")
     do i_atom=1,N_wz
        q_cell=q_cell + atom_wz(i_atom)%q
        p_cell(:)=p_cell(:) + atom_wz(i_atom)%rc(:)*atom_wz(i_atom)%q
!!$           format="(7X,i4,4x,a,3x,3(f17.12,3x),3x, f7.3 )"
!!$           write(output_unit,format) i_atom, atom_wz(i_atom)%name, atom_wz(i_atom)%rc(:),&
!!$                                             atom_wz(i_atom)%q
!!$           format="(a4,4x,3(f12.7,3x),5x,f10.7,3x,i4)"
!!$           write(scell_unit,format) atom_wz(i_atom)%name,atom_wz(i_atom)%rc(:), &
!!$                                    atom_wz(i_atom)%q, i_atom
     end do
!!$        format="(10x, 'Charge of symmetric cell :', f10.5 &
!!$               &/10x, 'Dipole of symmetric cell :', 3(f10.5,3x)/90(1h#))"
!!$        write(output_unit,format) q_cell,p_cell(:)
        !                             Make PC_array around the sym.cell
     call generate_lattice("direct",0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,R_max_pc,&
          print_vectors,full_print)

     N_ewald_pc=0
     do i_atom=1,N_wz
        lat: do i=1,lattice%num_vector
           rs(:)=atom_wz(i_atom)%rc(:)+lattice%v_cartesian(i,:)
           if(N_ewald_pc >= 1) then
              do j=1,N_ewald_pc
                 !if( all(rs==pc(j)%rc) ) then
                 if( all( abs(rs-pc(j)%rc)<=1.d-10) ) then
                    pc(j)%q=pc(j)%q + atom_wz(i_atom)%q
                    cycle lat
                 end if
              end do
           end if
           if(N_ewald_pc >= 1) then
              allocate(tmp_pc(N_ewald_pc), stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("MolMech: PC_ARRAY_MODULE: failed TMP_PC allocation")
              tmp_pc=pc
              deallocate(pc, stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("MolMech: PC_ARRAY_MODULE: failed PC deallocation")
              allocate(pc(N_ewald_pc+1),stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("MolMech: PC_ARRAY_MODULE: failed PC allocation")
              do k=1,N_ewald_pc
                 pc(k)=tmp_pc(k)
              end do
              deallocate(tmp_pc, stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("MolMech: PC_ARRAY_MODULE: failed TMP_PC deallocation")
           else
              allocate(pc(1),stat=alloc_stat)
              if(alloc_stat.ne.0) &
                   call error_handler("MolMech: PC_ARRAY_MODULE: failed PC(1) allocation")
           end if
           N_ewald_pc=N_ewald_pc+1
           pc(N_ewald_pc)%rc=rs
           pc(N_ewald_pc)%q=atom_wz(i_atom)%q
           pc(N_ewald_pc)%name=atom_wz(i_atom)%name
        end do lat
     end do
     !                                   Store Ewald_PC.xyz
!!$        write(output_unit,"(/37(1h#),' Ewald PC-array ',37(1h#)  )" )
!!$        write(output_unit,"(10X,'Ewald PC-array is stored in : ',a)") pc_xyz_file
!!$        write(pc_xyz_unit,"(i7/'comment line')" ) N_ewald_pc
!!$        write(output_unit,"(/30X,'X',17x,'Y',20x,'Z',17x,'Q' )")
     cluster(:)%V_pc=0.0_r8_kind
     q_tot=0.0_r8_kind
     do i=1,N_ewald_pc
        q_tot=q_tot+pc(i)%q
!!$           format="(7X,i4,4x,a,3x,3(f17.12,3x),3x, f7.3 )"
!!$           write(output_unit,format)  &
!!$                i, pc(i)%name, pc(i)%rc(:), pc(i)%q
!!$           format="(a4,4x,3(f12.7,3x),5x,f10.7,3x,i4)"
!!$           write(pc_xyz_unit,format) &
!!$                pc(i)%name,pc(i)%rc(:)+rc_screepAU,pc(i)%q, i
        do i_atom=1,n_cluster_atoms
           if( all(abs(pc(i)%rc-cluster(i_atom)%rc)<=1.d-1) ) cycle
           cluster(i_atom)%V_pc = cluster(i_atom)%V_pc + pc(i)%q / &
                sqrt(dot_product(pc(i)%rc-cluster(i_atom)%rc,pc(i)%rc-cluster(i_atom)%rc))
        end do
     end do

!!$        write(output_unit,"(/10x,'Total PC charge                  : ',f12.7)") q_tot
!!$        write(output_unit,"(/10X,'Coulomb potential on cluster atoms')")
!!$        format="(10x,i3,3x,a,4x,f12.7,4x,f12.7)"
!!$        do i_atom=1,n_cluster_atoms
!!$           write(output_unit,format) i_atom,cluster(i_atom)%name, cluster(i_atom)%V_pc
!!$        end do
!!!!!!!!!!!!!!!!!!!
     do i_atom=1,n_cluster_atoms
        cluster(i_atom)%V_pc=0.0_r8_kind
        k=0
        do i=1,N_ewald_pc
           R=sqrt(dot_product(pc(i)%rc-cluster(i_atom)%rc,pc(i)%rc-cluster(i_atom)%rc))
           if( R <=1.0e-1_r8_kind ) cycle
           k=k+1
!!$              write(output_unit,*) "R, i_atom, i, pc(i)%q / R  ",R, i_atom, i, pc(i)%q / R
              cluster(i_atom)%V_pc=cluster(i_atom)%V_pc + pc(i)%q / R
        end do
!!$           write(output_unit,*) "Cl  :", i_atom, cluster(i_atom)%V_pc, N_ewald_pc,k
     end do
!!!!!!!!!!!!!!!!!!

   end subroutine pc_generate
  !******************************************************************

  !******************************************************************
   subroutine Madelung(rl_in,V_madelung,n_centers)

     integer(i4_kind), intent(in)                       :: n_centers
     real(r8_kind),              dimension(n_centers,3) :: rl_in
     real(r8_kind), intent(out), dimension(n_centers)   :: V_madelung
     real(r8_kind),              dimension(3)           :: rc,v
     real(r8_kind)                                      :: zero=0.0_r8_kind
     real(r8_kind)                                      :: x_eps,R_max,G_max
     real(r8_kind)                                      :: R,G,Gr,M,eps=1.0e-15_r8_kind
     integer(i4_kind)                                   :: i_atom,i_vector,i_center,i

!!$        pi=acos(-1.d0)
     x_eps=r_f(myerfc,eps,zero)
     R_max=x_eps/eta                     ! upper estimate for "direct" lattice summation
     x_eps=r_f(expf,4*eta**2*eps,0.5_r8_kind)
     G_max=two*eta*x_eps                   ! upper estimate for "direct" lattice summation

!!$        eta=sqrt(((n_atoms_in_cell*pi**3)/V_cell**2)**(1d0/3d0))
!!$        x_eps=sqrt(-log(1.0d-8))
!!$        R_max=x_eps/eta
!!$        G_max=2*eta*x_eps

!!$        write(output_unit,*) 'Direct Sum R_max :', R_max,eta,x_eps
!!$        write(output_unit,*) 'Recipr.Sum G_max :', G_max, exp(-(G_max/(2*eta))**2)/G_max**2

     call generate_lattice("direct",zero,zero,zero,R_max,print_vectors)
!!$print*,lattice%num_vector
     do i_center=1,n_centers
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
              v(:)=lattice%v_cartesian(i_vector,:)+atomew(i_atom)%rc(:) - rc(:)
              R=sqrt(dot_product(v,v))
              if(R<=1.0e-15_r8_kind) R=zero
              if(R/=zero) then
                 M = M + myerfc(eta*R)/R
              else
                 M = M - 2*eta/sqrt(pi)
              end if
           end do
           V_madelung(i_center)=V_madelung(i_center) + atomew(i_atom)%q * M
        end do
     end do
     !                            Add the summ over reciprocal lattice

     call generate_lattice("reciprocal",zero,zero,zero,G_max,print_vectors,full_print)
!!$print*,lattice%num_vector
     do i_center=1,n_centers
        do i_atom=1,n_atoms_in_cell
           M=zero
           do i_vector=2,lattice%num_vector
              G=lattice%v_mod(i_vector)
              Gr=two*pi*dot_product(lattice%v_vector(i_vector,:), &
                   rl_in(i_center,:)-atomew(i_atom)%rl(:))
              M = M + exp(-(G/(two*eta))**2 )*cos(Gr) / G**2
           end do
           V_madelung(i_center)=V_madelung(i_center) + &
                atomew(i_atom)%q * (four*pi*M/V_cell - pi/(eta*eta*V_cell) )
        end do
     end do

   end subroutine Madelung
  !******************************************************************

  !******************************************************************
   function r_f(f,y,x0)  result(x)

     real(r8_kind), intent(in)       :: y,x0
     real(r8_kind)                   :: x,dx,x1,x2,tol=1.0e-15_r8_kind
     real(r8_kind), external         :: f

     dx=one
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
     x=(x1+x2)/two
     do while( dx > tol )
        if(f(x)> y) then
           x1=x
        else
           x2=x
        end if
        x=(x1+x2)/2.0_8
        dx=abs(x2-x1)
     end do
!!$     write(output_unit,*) 'eps  =', y
!!$     write(output_unit,*) 'F(X0)=',f(x),x
   end function r_f
   !******************************************************************

   !******************************************************************
   function expf(x)

     real(r8_kind)  :: expf,x

     expf=exp(-x**2)/x**2

   end function expf
   !******************************************************************

   !******************************************************************
   function myerfc(x)
     !     COMPLEMETARY ERROR FUNCTION USING CHEBYSHEV APPROXIMATION (YL LUKE
     !     1975 PP123-4) D,DD,SV AS IN PRESS, NUM. REC. 1988 "CHEBEV"
        integer(i4_kind), parameter :: na=25, nc=22
        real(r8_kind), dimension(0:na) :: a
        real(r8_kind), dimension(0:nc) :: c
        real(r8_kind)                  :: zero1=0.0, half1=0.5, one1=1.0,     &
             two1=2.0,three1=3.0, four1=4.0,     &
             sqpi2=1.128379167095513
        real(r8_kind)                  :: myerfc,x,z,d,dd,alpha,sv
        integer(i4_kind)               :: j
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
      d=zero1
      dd=zero1
      if (abs(x) < three1) then
         !         CALCULATE VIA ERF
         z=x/three1
         alpha=two1-four1*z*z
         do j=na,0,-1
            sv=d
            d=-alpha*d-dd+a(J)
            dd=sv
         end do
         myerfc=one1-sqpi2*z*(d-dd)
      else
         !      CALCULATE DIRECTLY
         z=abs(three1/x)
         alpha=two1-four1*z*z
         do j=nc,0,-1
            sv=d
            d=-alpha*d-dd+c(j)
            dd=sv
         end do
         if (x > zero1) then
            myerfc=half1*exp(-x*x)/x*(d+half1*alpha*dd)*sqpi2
         else
            myerfc=two1-half1*exp(-x*x)/(-x)*(d+half1*alpha*dd)*sqpi2
         end if
      end if

    end function myerfc
   !******************************************************************

   !******************************************************************
    subroutine inverse_matrix(A,n)

      integer(i4_kind), intent(in)                  :: n
      real(r8_kind), intent(inout),  dimension(n,n) :: a
      real(r8_kind) :: d,tol,amax,t,swap,&
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
!c        if(mod(NSWAP,2).ne.0) D=-D
      D=D*((-1.d0)**NSWAP)
      deallocate(ipv, stat=alloc_stat)
      if(alloc_stat /= 0) &
           call error_handler("Inverse_matrix [2]: Deallocation failed")

    END subroutine inverse_matrix
   !******************************************************************

   !******************************************************************
  end module pc_array_module






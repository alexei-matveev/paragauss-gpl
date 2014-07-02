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
!==================================================================
! Public interface of module
!==================================================================
module potential_calc_module
  !
  !  Purpose : Calculate a electrostatic potential of QM systems
  !            on choosen plane grid
  !
  !  Author: AS
  !  Date: 02/2001
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   03/2002
  ! Description: Calculation of Potential Derived Charges
  !              was implemented
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module, only: i4_kind, r8_kind
  use iounitadmin_module, only: output_unit, openget_iounit, &
       returnclose_iounit
  use potential_module, only: spacepoint_type, point_in_space, &
       N_points, V_pot_e, V_pot_n, V_pot_pc, dealloc_space_points, &
       deallocate_pot, start_read_poten_e, get_poten_n, bounds_free_poten, &
       send_recv_space_point, get_poten_pc
  use filename_module, only: inpfile
#ifdef FPP_AIX_XLF
  use matrix_module, only: matmult
# define MATMUL(a,b) matmult(a,b)
#endif
  implicit none
  save
  private

  !------------ public functions and subroutines ------------------
  public poten_calc_read, poten_calc_write, calc_plane_grid, grid2space_2d, &
       get_poten_and_shutdown_2d, calc_shell_grid, collect_poten_3d, &
       calc_poten_derive_charges
!==================================================================
! End of public interface of module
!==================================================================
  !== Interrupt end of public interface of module =================
  !-- Declaration of privat constant and variable  ----------------
  type poten
     real(kind=r8_kind) :: r(3)
     real(kind=r8_kind) :: pot
  end type poten

  type (poten), allocatable :: pot_on_plane(:,:)
  type (poten), allocatable :: pot_on_shell(:,:)

  logical, public :: esp_map, pdc

  character(len=10) :: potential_task
  logical, public :: use_saved_densmatrix
  logical, public :: V_electronic
  logical, public :: V_nuclear
  logical, public :: V_pc
  real(kind=r8_kind) :: first_point(3)
  real(kind=r8_kind) :: second_point(3)
  real(kind=r8_kind) :: third_point(3)
  real(kind=r8_kind) :: Xlimits(2)
  real(kind=r8_kind) :: Ylimits(2)
  integer(kind=i4_kind) :: Xgrid
  integer(kind=i4_kind) :: Ygrid
  real(kind=r8_kind) :: V_abs_limit
  character(len=4) :: output_units
  character(len=9) :: output_format

  character(len=10) :: df_potential_task= "pdc" ! "esp_map"
  logical :: df_use_saved_densmatrix = .false.
  logical :: df_V_electronic = .true.
  logical :: df_V_nuclear = .true.
  logical :: df_V_pc = .false.
  real(kind=r8_kind) :: df_first_point(3)= (/0.0,0.0,0.0/)
  real(kind=r8_kind) :: df_second_point(3)= (/1.0,0.0,0.0/)
  real(kind=r8_kind) :: df_third_point(3)= (/0.0,1.0,0.0/)
  real(kind=r8_kind) :: df_Xlimits(2)= (/0.0,1.0/)
  real(kind=r8_kind) :: df_Ylimits(2)= (/0.0,1.0/)
  integer(kind=i4_kind) :: df_Xgrid=50
  integer(kind=i4_kind) :: df_Ygrid=50
  real(kind=r8_kind) :: df_V_abs_limit=0.0_r8_kind
  character(len=4) :: df_output_units="a.u."
  character(len=9) :: df_output_format="gnuplot"

  integer(kind=i4_kind) :: n_shells
  real(kind=r8_kind) :: vdw_increments(8)
  character(len=6) :: pdc_constraint
  character(len=6) :: pdc_output
  integer(kind=i4_kind) :: n_fixed_atoms
  integer(kind=i4_kind) :: fixed_atoms(100)
  real(kind=r8_kind) :: fixed_charges(100)
  character(len=9) :: pdc_restraint
  real(kind=r8_kind) :: resp_a,resp_b
  integer(kind=i4_kind) :: n_symm_groups
  integer(kind=i4_kind) :: n_symm_atoms(30)
  integer(kind=i4_kind) :: symm_atoms(100)

  integer(kind=i4_kind) :: df_n_shells=4_i4_kind
  real(kind=r8_kind) :: df_vdw_increments(8)= (/1.0_r8_kind,1.2_r8_kind,1.4_r8_kind,1.6_r8_kind, &
       0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0.0_r8_kind/)
  character(len=6) :: df_pdc_constraint="charge" ! "dipole"
  character(len=6) :: df_pdc_output="short"
  integer(kind=i4_kind) :: df_n_fixed_atoms=0_i4_kind
  integer(kind=i4_kind) :: df_fixed_atoms(100)=0_i4_kind
  real(kind=r8_kind) :: df_fixed_charges(100)=0.0_r8_kind
  character(len=9) :: df_pdc_restraint="no" !"resp", "resp+symm"
  real(kind=r8_kind) :: df_resp_a=0.0005_r8_kind !0.001_r8_kind
  real(kind=r8_kind) :: df_resp_b=0.1_r8_kind
  integer(kind=i4_kind) :: df_n_symm_groups=0
  integer(kind=i4_kind) :: df_n_symm_atoms(30)=0
  integer(kind=i4_kind) :: df_symm_atoms(100)=0

  logical :: scilab, worksheet,gnuplot
  integer(kind=i4_kind) :: pdc_out_type
  logical,public :: charge_constr
  logical,public :: dipole_constr
  logical,public :: resp_restr
  logical,public :: symm_restr

  integer(kind=i4_kind) :: poligon_type,N_centers_on_sphere
  real(kind=r8_kind), allocatable :: r_sphere(:,:)

  real(kind=r8_kind), parameter :: au2ang = 0.529177249_r8_kind

  namelist/potential/ potential_task, use_saved_densmatrix, &  ! common options
       first_point,second_point,third_point, &             !
       Xlimits,Ylimits,Xgrid,Ygrid,V_electronic, &         ! ESP maps options
       V_nuclear,V_pc,V_abs_limit, &                       !
       output_units,output_format, &                       !
       n_shells,vdw_increments, &            !
       pdc_output,pdc_constraint, &          ! PDC options
       n_fixed_atoms, &                      !
       fixed_atoms,fixed_charges, &          !
       pdc_restraint, &                      !
       resp_a,resp_b, &                      !
       n_symm_groups, &                      !
       n_symm_atoms,symm_atoms               !
!------------ Subroutines -----------------------------------------

contains

!******************************************************************
  subroutine poten_calc_read()
    !------------ Modules used ------------------------------------
    use input_module
    use operations_module, only: operations_potential,operations_dipole
    !------------ Declaration of formal parameters ----------------
    !== End of interface ==========================================
    !------------ Declaration of local variables ------------------
    integer(kind=i4_kind) :: unit,status,i
    real(kind=r8_kind), parameter :: small=1.0e-03_r8_kind
    !------------ Executable code ---------------------------------

    potential_task=df_potential_task

    use_saved_densmatrix=df_use_saved_densmatrix

    V_electronic=df_V_electronic
    V_nuclear=df_V_nuclear
    first_point=df_first_point
    second_point=df_second_point
    third_point=df_third_point
    Xlimits=df_Xlimits
    Ylimits=df_Ylimits
    Xgrid=df_Xgrid
    Ygrid=df_Ygrid
    output_units=df_output_units
    output_format=df_output_format
    pdc_out_type=1

    n_shells=df_n_shells
    vdw_increments=df_vdw_increments
    n_fixed_atoms=df_n_fixed_atoms
    fixed_atoms=df_fixed_atoms
    fixed_charges=df_fixed_charges
    pdc_output=df_pdc_output
    pdc_constraint=df_pdc_constraint
    pdc_restraint=df_pdc_restraint
    resp_a=df_resp_a
    resp_b=df_resp_b
    n_symm_groups=df_n_symm_groups
    n_symm_atoms=df_n_symm_atoms
    symm_atoms=df_symm_atoms


    if(input_line_is_namelist("potential")) then
       unit = input_intermediate_unit()
       call input_read_to_intermediate
       read(unit, nml=potential, iostat=status)
       if (status .gt. 0) call input_error("calc_pot_read: namelist potential")
    endif

    if(operations_potential) then
       esp_map=(index(potential_task,"esp_map") /= 0 .or. index(potential_task,"ESP_MAP") /= 0)
       pdc=(index(potential_task,"pdc") /= 0 .or. index(potential_task,"PDC") /= 0)

       if(.not.esp_map .and. .not.pdc) &
            call error_handler("Potential_calc_module: Only esp_map or pdc tasks are possible to do")
       if(esp_map) then
          if(.not.V_electronic .and. .not.V_nuclear .and. .not. V_pc) &
               call error_handler("Potential_calc_module: Nothing to do !!!!")
          if(sqrt(dot_product(first_point-second_point,first_point-second_point)) <= small) &
               call error_handler("Potential_calc_module: first_point and second_point are coincided")
          if(sqrt(dot_product(first_point-third_point,first_point-third_point)) <= small) &
               call error_handler("Potential_calc_module: first_point and third_point are coincided")
          if(sqrt(dot_product(third_point-second_point,third_point-second_point)) <= small) &
               call error_handler("Potential_calc_module: third_point and second_point are coincided")
          if(abs(Xlimits(1)-Xlimits(2)) <= small) &
               call error_handler("Potential_calc_module: Limits on X axis are coinsided")
          if(abs(Ylimits(1)-Ylimits(2)) <= small) &
               call error_handler("Potential_calc_module: Limits on Y axis are coinsided")
          if(Xgrid <= 1) call error_handler("Potential_calc_module: Check xgrid value")
          if(Ygrid <= 1) call error_handler("Potential_calc_module: Check ygrid value")
          if(V_abs_limit < 0.0_r8_kind) call error_handler("Potential_calc_module: V_abs_limit < 0")
          if(index(output_units,"a.u.")==0 .and. index(output_units,"eV-a")==0) &
               call error_handler("Potential_calc_module: Only a.u. or eV-a units can be applied")
          if(index(output_format,"scilab")==0 .and. index(output_format,"SCILAB")==0 .and. &
               index(output_format,"worksheet")==0 .and. index(output_format,"WORKSHEET")==0 .and. &
               index(output_format,"gnuplot")==0 .and. index(output_format,"GNUPLOT")==0) &
               call error_handler( &
               "Potential_calc_module: Permited keywords are scilab(SCILAB), &
               & worksheet(WORKSHEET) or gnuplot(GNUPLOT)")
          scilab=.false.; worksheet=.false.; gnuplot=.false.
          if(index(output_format,"scilab")/=0 .or. index(output_format,"SCILAB")/=0) scilab=.true.
          if(index(output_format,"worksheet")/=0 .or. index(output_format,"WORKSHEET")/=0) worksheet=.true.
          if(index(output_format,"gnuplot")/=0 .or. index(output_format,"GNUPLOT")/=0) gnuplot=.true.
       else if(pdc) then
          charge_constr=.false.
          dipole_constr=.false.
          resp_restr=.false.
          symm_restr=.false.
          if(n_shells < 1) call error_handler("Potential_calc_module: n_shells has to be larger than 1")
          if(n_shells > 8) call error_handler("Potential_calc_module: n_shells has to be smaller than 8")
          if(index(pdc_output,"short") == 0 .and. index(pdc_output,"SHORT") == 0 .and. &
               index(pdc_output,"middle") == 0 .and. index(pdc_output,"MIDDLE") == 0 .and. &
               index(pdc_output,'long') == 0 .and. index(pdc_output,'LONG') == 0 ) &
               call error_handler( &
               "Potential_calc_module: pdc_output - only short(SHORT), middle(MIDDLE) &
               & or long(LONG) are possible")
          if(index(pdc_output,'short') /= 0 .or. index(pdc_output,'SHORT') /= 0) pdc_out_type=1
          if(index(pdc_output,'middle') /= 0 .or. index(pdc_output,'MIDDLE') /= 0) pdc_out_type=2
          if(index(pdc_output,'long') /= 0 .or. index(pdc_output,'LONG') /= 0 ) pdc_out_type=3
          charge_constr=(index(pdc_constraint,"charge") /= 0 .or. index(pdc_constraint,"CHARGE") /= 0)
          dipole_constr=(index(pdc_constraint,"dipole") /= 0 .or. index(pdc_constraint,"DIPOLE") /= 0)
          resp_restr=(index(pdc_restraint,"resp") /= 0 .or. index(pdc_restraint,"RESP") /= 0)
          if(resp_restr) then
             charge_constr=.true.
             dipole_constr=.false.
          end if
          symm_restr=(index(pdc_restraint,"resp+symm") /= 0 .or. index(pdc_restraint,"RESP+SYMM") /= 0)
          if(n_symm_groups==0) symm_restr=.false.
          if(symm_restr) then
             charge_constr=.true.
             resp_restr=.true.
             dipole_constr=.false.
          end if
          do i=1,n_symm_groups
             if(n_symm_atoms(i) <=1) n_symm_atoms(i)=0
          end do
          if(.not.symm_restr) n_symm_groups=0
          if(.not.charge_constr .and. .not.dipole_constr) &
               call error_handler( &
               "Potential_calc_module: pdc_constraint - only charge(CHARGE) &
               & or dipole(DIPOLE) are possible")
          if(dipole_constr) operations_dipole=.true.
       end if
    endif

  end subroutine poten_calc_read
!******************************************************************

!******************************************************************
  subroutine poten_calc_write(iounit)
    !------------ Modules used ------------------------------------
    use echo_input_module
    use operations_module, only: operations_potential,operations_echo_input_level
    !------------ Declaration of formal parameters ----------------
    integer(kind=i4_kind), intent(in) :: iounit
    !== End of interface ==========================================
    !------------ Declaration of local variables ------------------
    integer (i4_kind) :: i, nn
    !------------ Executable code ---------------------------------

    if( operations_potential ) then
       call start("POTENTIAL","POT_CALC_WRITE",iounit,operations_echo_input_level)
       call   word("POTENTIAL_TASK      ",potential_task      ,df_potential_task)
       call   flag("USE_SAVED_DENSMATRIX",use_saved_densmatrix,df_use_saved_densmatrix)
       if(esp_map) then
          call   flag("V_ELECTRONIC        ",v_electronic        ,df_v_electronic    )
          call   flag("V_NUCLEAR           ",v_nuclear           ,df_v_nuclear       )
          call   flag("V_PC                ",v_pc                ,df_v_pc            )
          call real_d("FIRST_POINT         ",first_point         ,df_first_point   ,3)
          call real_d("SECOND_POINT        ",second_point        ,df_second_point  ,3)
          call real_d("THIRD_POINT         ",third_point         ,df_third_point   ,3)
          call real_d("XLIMITS             ",xlimits             ,df_xlimits       ,2)
          call real_d("YLIMITS             ",ylimits             ,df_ylimits       ,2)
          call   intg("XGRID               ",xgrid               ,df_xgrid           )
          call   intg("YGRID               ",ygrid               ,df_ygrid           )
          call   real("V_ABS_LIMIT         ",v_abs_limit         ,df_v_abs_limit     )
          call   word("OUTPUT_UNITS        ",output_units        ,df_output_units    )
          call   word("OUTPUT_FORMAT       ",output_format       ,df_output_format   )
       else if(pdc) then
          call   intg("N_SHELLS            ",n_shells            ,df_n_shells        )
          call real_d("VDW_INCREMENTS      ",vdw_increments      ,df_vdw_increments,n_shells)
          call   intg("N_FIXED_ATOMS       ",n_fixed_atoms       ,df_n_fixed_atoms   )
          if(n_fixed_atoms > 0) then
             call intg_d("FIXED_ATOMS         ",fixed_atoms         ,df_fixed_atoms,n_fixed_atoms)
             call real_d("FIXED_CHARGES       ",fixed_charges       ,df_fixed_charges,n_fixed_atoms)
          end if
          call   word("PDC_CONSTRAINT      ",pdc_constraint      ,df_pdc_constraint  )
          call   word("PDC_RESTRAINT       ",pdc_restraint       ,df_pdc_restraint   )
          if(resp_restr) then
             call real("RESP_A              ",resp_a               ,df_resp_a)
             call real("RESP_B              ",resp_b               ,df_resp_b)
          end if
          if(symm_restr) then
             call   intg("N_SYMM_GROUPS       ",n_symm_groups  ,df_n_symm_groups)
             call intg_d("N_SYMM_ATOMS        ",n_symm_atoms   ,df_n_symm_atoms,n_symm_groups)
             nn=0
             do i=1,n_symm_groups
                nn=nn+n_symm_atoms(i)
             end do
             call intg_d("SYMM_ATOMS          ",symm_atoms     ,df_symm_atoms  ,nn )
          end if
          call   word("PDC_OUTPUT          ",pdc_output          ,df_pdc_output      )
       end if
       call stop()
    endif
!    contains
!    subroutine aaa
!    end subroutine aaa

  end subroutine poten_calc_write
!******************************************************************

!******************************************************************
  subroutine calc_plane_grid()
    !------------ Modules used ------------------------------------
    use unique_atom_module, only: N_unique_atoms,unique_atoms
    !------------ Declaration of formal parameters ----------------
    !== End of interface ==========================================
    !------------ Declaration of local variables ------------------
    real (r8_kind) :: center_coor(3), xaxis(3), yaxis(3), zaxis(3)
    real (r8_kind) :: length01, r01(3), rbuf(3), mat_rot(3,3)
    real(kind=r8_kind) :: atcoor(3,100),mat_rot_t(3,3)
    integer (i4_kind) :: iz(100), nat, atunit
    integer(kind=i4_kind) :: istat,i1,j1,k1
    !------------ Executable code ---------------------------------

    allocate(pot_on_plane(xgrid,ygrid),stat=istat)
    if(istat /= 0) call error_handler( &
         "potential_calc_module: allocate pot_on_plane failed")

    !saving X and Y points into ParaGauss INPUT directory
    if(scilab) call save_X_Y_points()

    !defining coordinate system on the grid plane
    length01=sqrt(dot_product(second_point-first_point,second_point-first_point))
    r01=second_point-first_point

    center_coor=first_point
    call vect_product(second_point-first_point,third_point-first_point,zaxis)
    zaxis=zaxis/sqrt(dot_product(zaxis,zaxis))

    xaxis=r01/sqrt(dot_product(r01,r01))

    call vect_product(zaxis,xaxis,yaxis)

    !build grid on choosen plane
    do i1=1,3
       mat_rot(i1,1)=xaxis(i1)/sqrt(dot_product(xaxis,xaxis))
       mat_rot(i1,2)=yaxis(i1)/sqrt(dot_product(yaxis,yaxis))
       mat_rot(i1,3)=zaxis(i1)/sqrt(dot_product(zaxis,zaxis))
    end do

    rbuf(3)=0.0_r8_kind
    do i1=1,xgrid
       rbuf(1)=Xlimits(1)+(i1-1)*(Xlimits(2)-Xlimits(1))/(xgrid-1)
       do j1=1,ygrid
          rbuf(2)=Ylimits(1)+(j1-1)*(Ylimits(2)-Ylimits(1))/(ygrid-1)
          pot_on_plane(i1,j1)%r=MATMUL(mat_rot,rbuf)+center_coor
       end do
    end do

    !converting atomic coordinates to the coordinate system connected
    !with the plane
    if(scilab) then
       mat_rot_t=transpose(mat_rot)
       k1=0
       do i1=1,N_unique_atoms
          do j1=1,unique_atoms(i1)%N_equal_atoms
             k1=k1+1
             rbuf=unique_atoms(i1)%position(:,j1)-center_coor
             atcoor(:,k1)=MATMUL(mat_rot_t,rbuf)
             iz(k1)=int(unique_atoms(i1)%Z)
          enddo
       enddo
       nat=k1

       atunit=openget_iounit(trim(inpfile('Cluster')), &
            form='formatted', status='unknown')
       write(atunit,'(i3)') nat
       do i1=1,nat
          if(trim(output_units) == "eV-a") atcoor(:,i1)=atcoor(:,i1)*au2ang
          write(atunit,'(i3,3(f17.8))') iz(i1),atcoor(:,i1)
       enddo
       call returnclose_iounit(atunit)
    endif

  contains
    subroutine vect_product(av,bv,cv)
      !------------ Declaration of formal parameters --------------
      real(kind=r8_kind),intent(in) :: av(3),bv(3)
      real(kind=r8_kind),intent(out) :: cv(3)
      !== End of interface ========================================
      !------------ Declaration of local variables ----------------
      !------------ Executable code -------------------------------
      cv(1)=av(2)*bv(3)-av(3)*bv(2)
      cv(2)=av(3)*bv(1)-av(1)*bv(3)
      cv(3)=av(1)*bv(2)-av(2)*bv(1)
    end subroutine vect_product
!------------------------------------------------------------------

!------------------------------------------------------------------
    subroutine save_X_Y_points()
      !------------ Declaration of formal parameters --------------
      !== End of interface ========================================
      !------------ Declaration of local variables ----------------
      real(kind=r8_kind) :: x,y
      integer(kind=i4_kind) :: grid_unit,xunit,yunit,units,i
      !------------ Executable code -------------------------------

      units=0
      if(index(output_units,"eV-a")/=0) units=1

      grid_unit=openget_iounit(trim(inpfile('xygrid')),  &
           form='formatted', status='unknown')
      write(grid_unit,'(i5)') xgrid
      write(grid_unit,'(i5)') ygrid
      write(grid_unit,'(i5)') units
      call returnclose_iounit(grid_unit)

      xunit=openget_iounit(trim(inpfile('Xpoints')),  &
           form='formatted', status='unknown')
      do i=1,xgrid
         x=Xlimits(1)+(i-1)*(Xlimits(2)-Xlimits(1))/(xgrid-1)
         if(index(output_units,"eV-a")/=0) x=x*au2ang
         write(xunit,'(f20.15)') x
      end do
      call returnclose_iounit(xunit)

      yunit=openget_iounit(trim(inpfile('Ypoints')),  &
           form='formatted', status='unknown')
      do i=1,ygrid
         y=Ylimits(1)+(i-1)*(Ylimits(2)-Ylimits(1))/(ygrid-1)
         if(trim(output_units) == "eV-a") y=y*au2ang
         write(yunit,'(f20.15)') y
      end do
      call returnclose_iounit(yunit)

    end subroutine save_X_Y_points

  end subroutine calc_plane_grid
!********************************************************************

!********************************************************************
  subroutine grid2space_2d()
    !
    ! Runs on all workers.
    !
    use comm, only: comm_rank
    use group_module, only: ylm_trafos, sub_group, group_coset, &
         symm_transformation_int, group_num_el, group_coset_decomp
    implicit none
    !== End of interface ============================================

    ! FIXME: why not everywhere?
    if (comm_rank () == 0) then
       call symm_sorted_grid ()
    endif

    call send_recv_space_point()

  contains
    subroutine symm_sorted_grid
      !------------ Declaration of formal parameters ----------------
      !== End of interface ==========================================
      !------------ Declaration of local variables ------------------
      integer(kind=i4_kind) :: n_equal,n_equal_check,N_total,n_size
      real(kind=r8_kind)    :: position(3),position2(3),position3(3)
      type(sub_group) :: local_groups
      type(group_coset) :: cosets
      type(symm_transformation_int) :: point_trafos
      logical, allocatable :: help_dim(:)
      type(spacepoint_type),allocatable :: buffer(:)

      real(kind=r8_kind),parameter :: small = 1.0e-10_r8_kind
      integer(kind=i4_kind) :: i,n,j,k,status,i1,j1,k1,l1

      !------------ Executable code ---------------------------------

      N_total=Xgrid*Ygrid

      allocate(help_dim(N_total),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: allocation HELP_DIM failed")
      help_dim=.false.

      allocate(buffer(N_total),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_modulr: allocation BUFFER failed")

      n_size=0
      i=0
      i1_xgrid: do i1=1,Xgrid
         j1_ygrid: do j1=1,Ygrid
            i=i+1
            help_d: if(.not.help_dim(i)) then
               n_size=n_size+1
               ! reorder coordinates of points as
               ! (x,y,z) --> (x,z,y) in order to comply with the
               ! convention for angular momentum l=1
               position(1) = pot_on_plane(i1,j1)%r(1)
               position(2) = pot_on_plane(i1,j1)%r(3)
               position(3) = pot_on_plane(i1,j1)%r(2)
               !
               ! determine local symmetry groups
               !
               ! now apply all symmetry operations to the position of the
               ! surface point
               n = 0
               do j=1,group_num_el
                  position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
                  if (dot_product(position2-position,position2-position) < small) then
                     n = n+1
                  endif
               enddo

               ! allocate group elements
               local_groups%num_el = n
               allocate(local_groups%elements(n))
               if ( status /= 0) call error_handler( &
                    "points_on_cavity_surface: allocation LOCAL_GROUPS  is failed")

               ! fill up group elements
               n = 0
               do j=1,group_num_el
                  position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
                  if (dot_product(position2-position,position2-position) < small) then
                     n = n+1
                     local_groups%elements(n) = j
                  end if
               enddo
               !
               ! now determine symmetry equivalent atoms
               !
               call group_coset_decomp(n_equal,local_groups,&
                    cosets,point_trafos%matrix)
               !
               ! search of positions of equal atoms
               !
               allocate(buffer(n_size)%position(3,n_equal), stat=status)
               if (status .ne. 0) call error_handler( &
                    "potential_calc_module: allocation of BUFFER%POSITION  failed")
               n_equal_check=0
               j_n_equal: do j=1,n_equal
                  position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets%elements(1,j)),position)
                  position3(1) = position2(1)
                  position3(2) = position2(3)
                  position3(3) = position2(2)
                  k=0
                  k1_xgrid: do k1=1,Xgrid
                     l1_ygrid: do l1=1,Ygrid
                        k=k+1
                        if(.not.help_dim(k)) then
                           if (sqrt(dot_product(position3-pot_on_plane(k1,l1)%r, &
                                position3-pot_on_plane(k1,l1)%r)) <= small) then

                              help_dim(k)=.true.
                              n_equal_check=n_equal_check+1
                              buffer(n_size)%position(:,j)=pot_on_plane(k1,l1)%r
                              buffer(n_size)%n_equal_points=n_equal
                           endif
                        endif
                     enddo l1_ygrid
                  enddo k1_xgrid
               enddo j_n_equal
               if (n_equal_check /= n_equal) call error_handler( &
                    "potential_calc_module: The program has calculated more equivalent points than &
                    & it had been able to find on the plane. Who fib from us? I do not know")

               deallocate(local_groups%elements, &
                    point_trafos%matrix,cosets%elements,stat=status)
               if (status .ne. 0 ) call error_handler( &
                    "potential_calc_module: deallocation of points HELPERS is failed")
            endif help_d
         enddo j1_ygrid
      enddo i1_xgrid

      N_points=n_size

      allocate(point_in_space(n_size),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: allocation POINT_IN_SPACE failed")

      do i=1,n_size
         point_in_space(i)%n_equal_points=buffer(i)%n_equal_points

         allocate(point_in_space(i)%position(3,buffer(i)%n_equal_points), stat=status)
         if (status .ne. 0 ) call error_handler( &
              "potential_calc_module: allocation of POINT_IN_SPACE(i)%position is failed")

         point_in_space(i)%position=buffer(i)%position

         deallocate(buffer(i)%position, stat=status)
         if (status .ne. 0 ) call error_handler( &
              "potential_calc_module: deallocation of BUFFER(i)%position is failed")
      enddo

      deallocate(buffer,stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: deallocation of BUFFER is failed")
      deallocate(help_dim,stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: deallocation of HELP_DIM  is failed")

    end subroutine symm_sorted_grid
  end subroutine grid2space_2d
!********************************************************************

!********************************************************************
  subroutine get_poten_and_shutdown_2d()
    !
    ! Executed on all workers.  Called from main_master(). Writes some
    ! output  to disk.   Historically this  procedure was  executed on
    ! master only  while slaves were spinning  in main_slave() waiting
    ! for orders.
    !
    use density_data_module, only: density_data_free
    use comm, only: comm_size, comm_rank
    !== End of interface ============================================

    real (r8_kind) :: x, y, conv
    integer (i4_kind) :: v_unit, i, j, k, l, m, istat

    if (V_electronic) then
       call start_read_poten_e()
    endif

    ! FIXME:  may  need some  tweaking  to  run  on all  workers.   In
    ! particular  wiriting  to  the  same  file  will  never  work  in
    ! parallel. Maybe just goto towards the end?
    if (comm_size() > 1) then
       ABORT("verify SPMD")
    endif
    ! FIXME: here or below?
    if (comm_rank() /= 0) goto 999 ! clean up and exit

    if (V_nuclear) call get_poten_n()
    if (V_pc) call get_poten_pc()

    if (index (potential_task, "esp_map") /= 0) then
       if (scilab) then
          v_unit = openget_iounit (trim (inpfile ('Vpoints')),  &
               form='formatted', status='unknown')
       else if (worksheet) then
          v_unit = openget_iounit (trim (inpfile ('V.txt')),  &
               form='formatted', status='unknown')
       else if (gnuplot) then
          v_unit = openget_iounit (trim (inpfile ('V.gpl')),  &
               form='formatted', status='unknown')
       endif
       conv = 1.0_r8_kind
       if (index (output_units, "eV-a") /= 0) conv = au2ang
       do i = 1, Xgrid
          if (worksheet .or. gnuplot) then
             x = (Xlimits(1) + (i - 1) * (Xlimits(2) - Xlimits(1)) / (xgrid - 1)) * conv
          endif
          do j = 1, Ygrid
             if (worksheet .or. gnuplot) then
                y = (Ylimits(1) + (j - 1) * (Ylimits(2) - Ylimits(1)) / (ygrid - 1)) * conv
             endif
             pot_on_plane(i, j) % pot = 0.0_r8_kind
             kk: do k = 1, N_points
                do l = 1, point_in_space(k) % N_equal_points
                   ! Maybe  norm2()  fortran  intrinsic  if  an  extra
                   ! sqrt() is not too expensive?:
                   associate (v => pot_on_plane(i, j) % r - point_in_space(k) % position(:, l))
                     if (dot_product (v, v) == 0.0) then
                        m = k
                        exit kk
                     endif
                   end associate
                enddo
             enddo kk
             if (V_electronic) then
                pot_on_plane(i, j) % pot = pot_on_plane(i, j) % pot + &
                     V_pot_e(m) / point_in_space(m) % N_equal_points
             endif
             if (V_nuclear) then
                pot_on_plane(i, j) % pot= pot_on_plane(i, j) % pot + &
                     V_pot_n(m) / point_in_space(m) % N_equal_points
             endif
             if (V_pc) then
                pot_on_plane(i, j) % pot = pot_on_plane(i, j) % pot + &
                     V_pot_pc(m) / point_in_space(m) % N_equal_points
             endif
             if (V_abs_limit > 0.0_r8_kind .and. &
                  abs (pot_on_plane(i, j) % pot) > V_abs_limit) then
                pot_on_plane(i, j) % pot = sign (V_abs_limit, pot_on_plane(i, j) % pot)
             end if
             if (index (output_units, "eV-a") /= 0) then
                pot_on_plane(i, j) % pot = pot_on_plane(i, j) % pot * 27.211652_r8_kind
             endif
             if (scilab) then
                write (v_unit, '(f25.15)') pot_on_plane(i, j) % pot
             else if (worksheet .or. gnuplot) then
                write (v_unit, '(3f25.15)') x, y, pot_on_plane(i, j) % pot
                if (gnuplot .and. j == Ygrid) write (v_unit, '(a)') ""
             endif
          enddo
       enddo
       call returnclose_iounit (v_unit)

       deallocate (pot_on_plane, stat=istat)
       if ( istat /= 0) call error_handler( &
              "potential_calc_module: deallocation pot_on_plane failed")
    endif

    ! FIXME: I assume this staff should be executed on all workers?
999 continue
    call dealloc_space_points() ! no comm
    call deallocate_pot()       ! no comm

    if (V_electronic) then
       call bounds_free_poten() ! no comm
       call density_data_free() ! no comm
    endif
  end subroutine get_poten_and_shutdown_2d
!********************************************************************

!********************************************************************
  subroutine calc_shell_grid()
    !
    ! Runs on all workers.
    !
    use comm, only: comm_rank
    use atoms_data_module
    use group_module, only: sub_group, group_coset, &
         symm_transformation_int, group_coset_decomp
    use unique_atom_module, only: N_unique_atoms,unique_atoms
    use symmetry_data_module, only: symmetry_data_point_group
    use symm_module, only: symm_adapt_centers
    implicit none
    !== End of interface ============================================

    integer (i4_kind) :: N_spheres, n_rotations
    integer(kind=i4_kind),allocatable :: n_points_on_shell(:)
    real(kind=r8_kind) :: radius
    real(kind=r8_kind),allocatable :: xyz_sphere(:,:)
    real (r8_kind), allocatable :: xyz_vertex(:, :)
    integer (i4_kind), allocatable :: tri_index(:,:)
    real(kind=r8_kind),allocatable :: xyz_centers(:,:)
    character(len=4) :: name_point_group

    ! FIXME: why not let everyone do the same?
    if (comm_rank () /= 0) goto 999 ! receive results

    call calc_shells()          ! no comm

    name_point_group =symmetry_data_point_group()

    if ( (name_point_group=='C1  ')  .or. &
         (name_point_group=='C2  ')  .or. &
         (name_point_group=='C3  ')  .or. &
         (name_point_group=='C5  ')  .or. &
         (name_point_group=='CS  ')  .or. &
         (name_point_group=='Ci  ')  .or. &
         (name_point_group=='C2V ')  .or. &
         (name_point_group=='C2H ')  .or. &
         (name_point_group=='C3V ')  .or. &
         (name_point_group=='C5V ')  .or. &
         (name_point_group=='S6  ')  .or. &
         (name_point_group=='S10 ')  .or. &
         (name_point_group=='D2  ')  .or. &
         (name_point_group=='D2H ')  .or. &
         (name_point_group=='D3  ')  .or. &
         (name_point_group=='D3D ')  .or. &
         (name_point_group=='D5  ')  .or. &
         (name_point_group=='D5D ')  .or. &
         (name_point_group=='I   ')  .or. &
         (name_point_group=='IH  ')  ) then
       call generate_dodecahedron()
       poligon_type=1
    elseif ( (name_point_group=='C4  ') .or. &
         (name_point_group=='C4V ')  .or. &
         (name_point_group=='C4H ')  .or. &
         (name_point_group=='D4  ')  .or. &
         (name_point_group=='D4H ')  .or. &
         (name_point_group=='S4  ')  .or. &
         (name_point_group=='D2D ')  .or. &
         (name_point_group=='O   ')  .or. &
         (name_point_group=='OH  ')  .or. &
         (name_point_group=='T   ')  .or. &
         (name_point_group=='TH  ')  .or. &
         (name_point_group=='TD  ') ) then
       call generate_cube()
        poligon_type=2
   else
       n_rotations=0
       if( (name_point_group=='C3H ')  .or. &
           (name_point_group=='D3H ') ) then
          n_rotations=3
       elseif( (name_point_group=='C5H ')  .or. &
            (name_point_group=='D5H ') ) then
          n_rotations=5
       elseif( (name_point_group=='C6  ')  .or. &
            (name_point_group=='C6V ')  .or. &
            (name_point_group=='C6H ')  .or. &
            (name_point_group=='D6  ')  .or. &
            (name_point_group=='D6H ') ) then
          n_rotations=6
       elseif( (name_point_group=='C7  ')  .or. &
            (name_point_group=='C7V ')  .or. &
            (name_point_group=='C7H ')  .or. &
            (name_point_group=='D7  ')  .or. &
            (name_point_group=='D7H ') ) then
          n_rotations=7
       elseif( (name_point_group=='C8  ')  .or. &
            (name_point_group=='C8V ')  .or. &
            (name_point_group=='C8H ')  .or. &
            (name_point_group=='D4D ')  .or. &
            (name_point_group=='S8  ')  .or. &
            (name_point_group=='D8  ')  .or. &
            (name_point_group=='D8H ') ) then
          n_rotations=8
       elseif( (name_point_group=='C9  ')  .or. &
            (name_point_group=='C9V ')  .or. &
            (name_point_group=='C9H ')  .or. &
            (name_point_group=='D9  ')  .or. &
            (name_point_group=='D9H ') ) then
          n_rotations=9
       elseif( (name_point_group=='C10 ')  .or. &
            (name_point_group=='C10V')  .or. &
            (name_point_group=='C10H')  .or. &
            (name_point_group=='D10 ')  .or. &
            (name_point_group=='D10H') ) then
          n_rotations=10
       elseif( (name_point_group=='S12 ')  .or. &
            (name_point_group=='D6D ')) then
          n_rotations=12
       elseif( (name_point_group=='S14 ')  .or. &
            (name_point_group=='D7D ')) then
          n_rotations=14
       elseif( (name_point_group=='S16 ')  .or. &
            (name_point_group=='D8D ')) then
          n_rotations=16
       elseif( (name_point_group=='S18 ')  .or. &
            (name_point_group=='D9D ')) then
          n_rotations=18
       elseif( (name_point_group=='S20 ')  .or. &
            (name_point_group=='D10D')) then
          n_rotations=20
       endif
       if(n_rotations==0) then
          call error_handler( &
               "potential_calc_module: S12 ... S20 and  &
               & D6D ... D10D are yet missing &
               & still .Sorry")
       endif
       call generate_doublepyramide(n_rotations)
       poligon_type=3
    endif

    call put_points_on_shells() ! no comm

    call symm_sorted_grid()     ! no comm

999 continue
    call send_recv_space_point() ! does comm

  contains
    !------------------------------------------------------------
    subroutine calc_shells()
      real(kind=r8_kind) :: v_d_w_r
      integer(kind=i4_kind) :: vdW_index
      integer(kind=i4_kind) :: i,j,k,l,status
!!$      integer(i4_kind) :: m !@@@@@@@@@@@

      N_spheres = 0
      do i=1,N_unique_atoms
         N_spheres = N_spheres + unique_atoms(i)%N_equal_atoms
      enddo
!!$      do i=1,pointcharge_N                                             !@@@@@@@@@
!!$         N_spheres = N_spheres + pointcharge_array(i)%N_equal_charges  !@@@@@@@@@
!!$      end do                                                           !@@@@@@@@@

      allocate(r_sphere(N_spheres,n_shells), stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: allocation of r_sphere is failed")
      allocate(xyz_sphere(N_spheres,3), stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: allocation of xyz_sphere is failed")

      k=0
      n_uni: do i=1,N_unique_atoms
         vdW_index=int(unique_atoms(i)%Z)
         v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind
         do j=1,unique_atoms(i)%N_equal_atoms
            k=k+1
            do l=1,n_shells
               r_sphere(k,l)=v_d_w_r*vdw_increments(l)/au2ang
            enddo
            xyz_sphere(k,:)=unique_atoms(i)%position(:,j)
         enddo
      enddo n_uni
!!$      n_pc: do i=1,pointcharge_N                                        !@@@@@@@@@
!!$         do m=1,98                                                      !@@@@@@@@@
!!$            if(index(pointcharge_array(i)%name,atom_name(m)) /= 0) exit !@@@@@@@@@
!!$         end do                                                         !@@@@@@@@@
!!$         vdW_index=m                                                    !@@@@@@@@@
!!$         v_d_w_r=R_def_rap(vdW_index)/2.0_r8_kind                       !@@@@@@@@@
!!$         do j=1,pointcharge_array(i)%N_equal_charges                    !@@@@@@@@@
!!$            k=k+1                                                       !@@@@@@@@@
!!$            do l=1,n_shells                                             !@@@@@@@@@
!!$               r_sphere(k,l)=v_d_w_r*vdw_increments(l)/au2ang           !@@@@@@@@@
!!$            enddo                                                       !@@@@@@@@@
!!$            xyz_sphere(k,:)=pointcharge_array(i)%position(:,j)          !@@@@@@@@@
!!$         enddo                                                          !@@@@@@@@@
!!$      end do n_pc                                                       !@@@@@@@@@

    end subroutine calc_shells
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine generate_cube()
      !generate set of surface points on a sphere starting from a cube
      !and do subdivisions
      !** End of interface *****************************************
      real(kind=r8_kind) :: a,d
      real(kind=r8_kind) :: r_buf,xyz_buf(3)
      integer(kind=i4_kind) :: tang_buf(4),local_point_factor
      integer(kind=i4_kind) :: N_dim_ver,N_dim_ver_st,N_dim_tri,N_dim_tri_st
      integer(kind=i4_kind) :: N_dim_ver_next,N_dim_tri_next
      integer(kind=i4_kind) :: i,j,k,l,m,n,l1,n1,ind,i1,status

      a=sqrt(2.0_r8_kind)/2.0_r8_kind
      radius=sqrt(1.5_r8_kind)
      d=(radius-a)**2+1.0_r8_kind

      N_dim_tri=24
      N_dim_tri_st=24
      N_dim_ver=14
      N_dim_ver_st=14

      local_point_factor=1 !!!!!!!!!!!

      i=1
      do
       if(i>local_point_factor) exit
         N_dim_ver=N_dim_ver+N_dim_tri*3/2
         N_dim_tri=N_dim_tri*4
         if(i==1) then
            N_dim_ver_next=N_dim_ver
            N_dim_tri_next=N_dim_tri
         endif
         i=i+1
      enddo

      allocate(xyz_vertex(N_dim_ver,3),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_cube: allocation of xyz_vertex is failed")


      !definition of initial set of points
      xyz_vertex=0.0_r8_kind

      xyz_vertex(1:4,1)=a
      xyz_vertex(5:8,1)=-a
      xyz_vertex(9,1)=radius
      xyz_vertex(10,1)=-radius
      xyz_vertex(11:14,1)=0.0_r8_kind

      xyz_vertex(1,2)=-a
      xyz_vertex(2,2)=a
      xyz_vertex(3,2)=-a
      xyz_vertex(4,2)=a
      xyz_vertex(5,2)=-a
      xyz_vertex(6,2)=a
      xyz_vertex(7,2)=-a
      xyz_vertex(8,2)=a
      xyz_vertex(9:10,2)=0.0_r8_kind
      xyz_vertex(11,2)=radius
      xyz_vertex(12,2)=-radius
      xyz_vertex(13:14,2)=0.0_r8_kind

      xyz_vertex(1:2,3)=a
      xyz_vertex(3:4,3)=-a
      xyz_vertex(5:6,3)=a
      xyz_vertex(7:8,3)=-a
      xyz_vertex(9:12,3)=0.0_r8_kind
      xyz_vertex(13,3)=radius
      xyz_vertex(14,3)=-radius

      !definition of triangles

      allocate(tri_index(N_dim_tri,3),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_cube: allocation of tri_index is failed")
      tri_index=0

      ind=1
      do i=1,6
         k=i+8
         m=4*(i-1)+1
         tri_index(m:m+3,1)=k

         n=1
         do j=1,8
            xyz_buf=xyz_vertex(k,:)-xyz_vertex(j,:)
            r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
            if (r_buf <= sqrt(d)+0.01_r8_kind) then
               tang_buf(n)=j
               n=n+1
               if (n > 4) exit
            endif

         enddo

         do l=1,3
            l1=tang_buf(l)
            do n=l+1,4
               n1=tang_buf(n)
               xyz_buf=xyz_vertex(n1,:)-xyz_vertex(l1,:)
               r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
               if (r_buf <= 2*a+0.01_r8_kind) then
                  tri_index(ind,2)=l1
                  tri_index(ind,3)=n1
                  ind=ind+1
               endif
            enddo
         enddo
      enddo

      !dividing each triangles in four new triangles
      do j=1,local_point_factor
         call more_triangles(N_dim_tri_next,N_dim_ver_st,N_dim_tri_st)
         N_dim_ver_st=N_dim_ver_next
         N_dim_tri_st=N_dim_tri_next
         N_dim_ver_next= N_dim_ver_st+N_dim_tri_st*3/2
         N_dim_tri_next= N_dim_tri_st*4
      enddo

      ! definition of triangle centers
      allocate(xyz_centers(N_dim_tri,3),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_cube: allocation of xyz_centers is failed")

      do i=1,N_dim_tri_st
         i1=tri_index(i,1)
         l1=tri_index(i,2)
         n1=tri_index(i,3)
         xyz_buf=(xyz_vertex(i1,:)+xyz_vertex(l1,:)+xyz_vertex(n1,:))/3.0_r8_kind
         xyz_centers(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
      enddo

      N_centers_on_sphere=N_dim_tri_st

      deallocate(xyz_vertex,tri_index,stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_cube: deallocation of xyz_vertex is failed")

    end subroutine generate_cube
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine generate_dodecahedron()
      !generate points due to a pentakis-dodecahedron
      !rotate in right symmetry
      !subdivide if needed
      !** End of interface *****************************************
      real(kind=r8_kind) :: x,y,z
      real(kind=r8_kind) :: t,s,h,r,b,c1,c2,h1,h2,rr,cosT,sinT
      real(kind=r8_kind) :: r_buf,xyz_buf(3)
      integer(kind=i4_kind) :: tang_buf(5),local_point_factor
      integer(kind=i4_kind) :: N_dim_ver,N_dim_tri,N_dim_ver_st,N_dim_tri_st
      integer(kind=i4_kind) :: N_dim_ver_next,N_dim_tri_next
      real(kind=r8_kind), parameter :: cos12=0.8660254037844_r8_kind
      real(kind=r8_kind), parameter :: sin12=0.5_r8_kind
      real(kind=r8_kind), parameter :: cos20=0.9510565162952_r8_kind
      real(kind=r8_kind), parameter :: sin20=0.3090169943749_r8_kind
      integer (i4_kind) :: i, j, k, l, m, n, l1, n1, ind, i1, status

      !definition a pentakis dodecahedron
      t=(sqrt(5.0_r8_kind)+1.0_r8_kind)/2.0_r8_kind
      s=sqrt(3.0_r8_kind-t)/2.0_r8_kind
      h=(2.0_r8_kind*s+(t+1.0_r8_kind)*sqrt(t+2.0_r8_kind))/4.0_r8_kind
      radius=sqrt(h**2+s**2)
      b=radius*sqrt(3.0_r8_kind)/3.0_r8_kind
      c1=(2.0_r8_kind*b+h)/5.0_r8_kind
      c2=(2.0_r8_kind*h+2.0_r8_kind*b+s)/5.0_r8_kind
      rr=sqrt(c1**2+c2**2)
      h1=c1*radius/rr
      h2=c2*radius/rr
      cosT=h2/sqrt(h1**2+h2**2)
      sinT=sqrt(1.0_r8_kind-cosT**2)
      r=sqrt(radius**2-1.0_r8_kind)

      N_dim_tri=60
      N_dim_tri_st=60
      N_dim_ver=32
      N_dim_ver_st=32

      local_point_factor=1 !!!!!!!!!!!!!!

      i=1
      do
       if(i>local_point_factor-1) exit
         N_dim_ver=N_dim_ver+(N_dim_tri*3)/2
         N_dim_tri=N_dim_tri*4
         if(i==1) then
            N_dim_ver_next=N_dim_ver
            N_dim_tri_next=N_dim_tri
         endif
         i=i+1
      enddo

      allocate(xyz_vertex(N_dim_ver,3),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_dodecahedron: allocation of xyz_vertex is failed")

      !definition of initial set of points
      xyz_vertex=0.0_r8_kind

      xyz_vertex(1:4,1)=0.0_r8_kind
      xyz_vertex(5,1)=s
      xyz_vertex(7,1)=s
      xyz_vertex(6,1)=-s
      xyz_vertex(8,1)=-s
      xyz_vertex(9:10,1)=h
      xyz_vertex(11:12,1)=-h
      xyz_vertex(13,1)=b
      xyz_vertex(15,1)=b
      xyz_vertex(17,1)=b
      xyz_vertex(19,1)=b
      xyz_vertex(14,1)=-b
      xyz_vertex(16,1)=-b
      xyz_vertex(18,1)=-b
      xyz_vertex(20,1)=-b
      xyz_vertex(21,1)=h1
      xyz_vertex(22,1)=-h1
      xyz_vertex(23,1)=h1
      xyz_vertex(24,1)=-h1
      xyz_vertex(25:29,1)=0.0_r8_kind
      xyz_vertex(29:30,1)=h2
      xyz_vertex(31:32,1)=-h2

      xyz_vertex(1,2)=-s
      xyz_vertex(3,2)=-s
      xyz_vertex(2,2)=s
      xyz_vertex(4,2)=s
      xyz_vertex(5:6,2)=h
      xyz_vertex(7:8,2)=-h
      xyz_vertex(9:12,2)=0.0_r8_kind
      xyz_vertex(13:14,2)=b
      xyz_vertex(17:18,2)=b
      xyz_vertex(15:16,2)=-b
      xyz_vertex(19:20,2)=-b
      xyz_vertex(21:24,2)=0.0_r8_kind
      xyz_vertex(25,2)=-h2
      xyz_vertex(26:27,2)=h2
      xyz_vertex(28,2)=-h2
      xyz_vertex(29,2)=-h1
      xyz_vertex(30,2)=h1
      xyz_vertex(31,2)=-h1
      xyz_vertex(32,2)=h1

      xyz_vertex(1:2,3)=h
      xyz_vertex(3:4,3)=-h
      xyz_vertex(5:8,3)=0.0_r8_kind
      xyz_vertex(9,3)=s
      xyz_vertex(11,3)=s
      xyz_vertex(10,3)=-s
      xyz_vertex(12,3)=-s
      xyz_vertex(13:16,3)=b
      xyz_vertex(17:20,3)=-b
      xyz_vertex(21:22,3)=h2
      xyz_vertex(23:24,3)=-h2
      xyz_vertex(25:26,3)=h1
      xyz_vertex(27:28,3)=-h1
      xyz_vertex(29:32,3)=0.0_r8_kind

      !definition of triangles

      allocate (tri_index(N_dim_tri,3), stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_dodecahedron: allocation of tri_index is failed")
      tri_index = 0

      ind=1
      do i = 1, 12
         k = i + 20
         m=5*(i-1)+1
         tri_index(m:m+4,1) = k

         n=1
         do j=1,20
            xyz_buf=xyz_vertex(k,:)-xyz_vertex(j,:)
            r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
            if (r_buf <= 1.06_r8_kind) then
               tang_buf(n)=j
               n=n+1
               if (n > 5) exit
            endif

         enddo

         do l=1,4
            l1 = tang_buf(l)
            do n=l+1,5
               n1 = tang_buf(n)
               xyz_buf=xyz_vertex(n1,:)-xyz_vertex(l1,:)
               r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
               if (r_buf <= 2.0_r8_kind*s+0.001_r8_kind) then
                  tri_index(ind, 2) = l1
                  tri_index(ind, 3) = n1
                  ind=ind+1
               endif
            enddo
         enddo
      enddo

      ! orientation of points respect to symmetry point group
      if ( (name_point_group=='C1  ') .or. &
           (name_point_group=='C2  ') .or. &
           (name_point_group=='CS  ') .or. &
           (name_point_group=='C2V ') .or. &
           (name_point_group=='C2H ') .or. &
           (name_point_group=='I   ') .or. &
           (name_point_group=='IH  ') .or. &
           (name_point_group=='CI  ')  .or. &
           (name_point_group=='D2  ')  .or. &
           (name_point_group=='D2H ') ) then
      endif
      if ( (name_point_group=='C3  ') .or. &
           (name_point_group=='C3V ') .or. &
           (name_point_group=='D3  ')  .or. &
           (name_point_group=='D3D ')  .or. &
           (name_point_group=='S6  ')  ) then
         do i=1,N_dim_ver_st
            x=xyz_vertex(i,1)*s/radius+xyz_vertex(i,3)*h/radius
            z=-xyz_vertex(i,1)*h/radius+xyz_vertex(i,3)*s/radius
            xyz_vertex(i,1)=x
            xyz_vertex(i,3)=z
         enddo
         if ((name_point_group=='D3  ')  .or. &
             (name_point_group=='D3D ') ) then
            do i=1,N_dim_ver_st
               x=cos12*xyz_vertex(i,1)-sin12*xyz_vertex(i,2)
               y=sin12*xyz_vertex(i,1)+cos12*xyz_vertex(i,2)
               xyz_vertex(i,1)=x
               xyz_vertex(i,2)=y
            enddo
         end if
      endif
      if ( (name_point_group=='C5  ') .or. &
           (name_point_group=='C5V ') .or. &
           (name_point_group=='S10 ')  .or. &
           (name_point_group=='D5  ')  .or. &
           (name_point_group=='D5D ')  .or. &
           (name_point_group=='D5H ') ) then
         do i=1,N_dim_ver_st
            x=xyz_vertex(i,1)*cosT+xyz_vertex(i,3)*sinT
            z=-xyz_vertex(i,1)*sinT+xyz_vertex(i,3)*cosT
            xyz_vertex(i,1)=x
            xyz_vertex(i,3)=z
         enddo
         if ((name_point_group=='D5  ') .or. &
             (name_point_group=='D5D '))  then
            do i=1,N_dim_ver_st
               x=cos20*xyz_vertex(i,1)-sin20*xyz_vertex(i,2)
               y=sin20*xyz_vertex(i,1)+cos20*xyz_vertex(i,2)
               xyz_vertex(i,1)=x
               xyz_vertex(i,2)=y
            enddo
         endif
      endif

      do j=1,local_point_factor-1
         call more_triangles(N_dim_tri_next,N_dim_ver_st,N_dim_tri_st)
         N_dim_ver_st=N_dim_ver_next
         N_dim_tri_st=N_dim_tri_next
         N_dim_ver_next= N_dim_ver_st+(N_dim_tri_st*3)/2
         N_dim_tri_next= N_dim_tri_st*4
      enddo

      ! definition of triangle centers
      allocate(xyz_centers(N_dim_tri,3),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_dodecahedron: allocation of xyz_centers is failed")

      do i=1,N_dim_tri_st
         i1 = tri_index(i, 1)
         l1 = tri_index(i, 2)
         n1 = tri_index(i, 3)
         xyz_buf=(xyz_vertex(i1,:)+xyz_vertex(l1,:)+xyz_vertex(n1,:))/3.0_r8_kind
         xyz_centers(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
      enddo

      N_centers_on_sphere=N_dim_tri_st

      deallocate (xyz_vertex, tri_index, stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_dodecahedron: deallocation of xyz_vertex is failed")

    end subroutine generate_dodecahedron
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine generate_doublepyramide(n_rotations)
      ! simple extension for all other point groups,
      ! not very well suited, because surface areas my be
      ! very acute angled triangles
      !** End of interface *****************************************
      real(kind=r8_kind) :: xyz_buf(3)
      integer(kind=i4_kind) :: N_dim_ver,N_dim_tri,N_dim_ver_st,N_dim_tri_st
      integer(kind=i4_kind) :: N_dim_ver_next,N_dim_tri_next
      integer(kind=i4_kind) :: i,j,l1,n1,i1,status
      integer(kind=i4_kind) :: n_rotations,local_point_factor
      real(kind=r8_kind) :: alph, r_n_rot
      real(kind=r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind
      real(kind=r8_kind), parameter :: cos12=0.8660254037844_r8_kind
      real(kind=r8_kind), parameter :: sin12=0.5_r8_kind

      radius=1.0_r8_kind

      N_dim_tri=2*n_rotations
      N_dim_tri_st=2*n_rotations
      N_dim_ver=n_rotations+2
      N_dim_ver_st=n_rotations+2

      local_point_factor=2

      i=1
      do
       if(i>local_point_factor) exit
         N_dim_ver=N_dim_ver+N_dim_tri*3/2
         N_dim_tri=N_dim_tri*4
         if(i==1) then
            N_dim_ver_next=N_dim_ver
            N_dim_tri_next=N_dim_tri
         endif
         i=i+1
      enddo

      allocate(xyz_vertex(N_dim_ver,3),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_doublepyramide: allocation of xyz_vertex is failed")

      !definition of initial set of points
      xyz_vertex=0.0_r8_kind

      xyz_vertex(1,1:2)=0.0_r8_kind
      xyz_vertex(1,3)=radius
      xyz_vertex(N_dim_ver_st,1:2)=0.0_r8_kind
      xyz_vertex(N_dim_ver_st,3)=-radius
      xyz_vertex(2:N_dim_ver_st-1,3)=0.0_r8_kind

      if (n_rotations == 3) then
        xyz_vertex(2,1)=1.0_r8_kind
        xyz_vertex(2,2)=0.0_r8_kind
        xyz_vertex(3,1)=-sin12
        xyz_vertex(3,2)=cos12
        xyz_vertex(4,1)=-sin12
        xyz_vertex(4,2)=-cos12
      else
        r_n_rot=real(n_rotations,kind=r8_kind)
        do i=0,n_rotations-1
          alph = (real(i,kind=r8_kind))/r_n_rot*2.0_r8_kind*pi
          xyz_vertex(i+2,1)=cos(alph)
          xyz_vertex(i+2,2)=sin(alph)
        enddo
      endif

      allocate (tri_index(N_dim_tri, 3), stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_doublepyramide: allocation of tri_index is failed")
      tri_index = 0

      do i = 1, n_rotations
        tri_index(i, 1) = i + 1
        tri_index(i, 2) = i + 2
        tri_index(i,3) = 1
        tri_index(n_rotations + i, 1) = i + 1
        tri_index(n_rotations + i, 2) = i + 2
        tri_index(n_rotations + i, 3) = N_dim_ver_st
      enddo
      tri_index(n_rotations, 2) = 2
      tri_index(N_dim_tri_st, 2) = 2


      !dividing each triangles in four new triangles
      do j=1,local_point_factor
         call more_triangles(N_dim_tri_next,N_dim_ver_st,N_dim_tri_st)
         N_dim_ver_st=N_dim_ver_next
         N_dim_tri_st=N_dim_tri_next
         N_dim_ver_next= N_dim_ver_st+N_dim_tri_st*3/2
         N_dim_tri_next= N_dim_tri_st*4
      enddo


      ! definition of triangle centers
      allocate(xyz_centers(N_dim_tri,3),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_doublepyramide: allocation of xyz_centers is failed")

      do i = 1, N_dim_tri_st
         i1 = tri_index(i, 1)
         l1 = tri_index(i, 2)
         n1 = tri_index(i, 3)
         xyz_buf=(xyz_vertex(i1,:)+xyz_vertex(l1,:)+xyz_vertex(n1,:))/3.0_r8_kind
         xyz_centers(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
      enddo

      N_centers_on_sphere=N_dim_tri_st

      deallocate (xyz_vertex, tri_index, stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:generate_doublepyramide: deallocation of xyz_vertex is failed")

    end subroutine generate_doublepyramide
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine more_triangles(n_ind,n_xyz_st,n_ind_st)
      !subdivision routine, each triangle is devided in four by
      !halving the edges
      !** End of interface *****************************************

      integer(kind=i4_kind), intent(in) :: n_ind
      integer(kind=i4_kind), intent(inout) :: n_xyz_st,n_ind_st

      real(kind=r8_kind) :: xyz_t_buf(3),xyz_tt(3)
      integer (i4_kind), allocatable :: ind_buf(:,:)
      integer(kind=i4_kind) :: new_numbers(6)
      integer(kind=i4_kind) :: status,next_i,neighbour,i1,i2
      integer(kind=i4_kind) :: i,j,k,l,m,m1,n,n1
      real(kind=r8_kind), parameter :: small=1.0e-11_r8_kind

      allocate(ind_buf(n_ind,3),stat=status)
      if ( status /= 0) call error_handler( &
           "more_triangles: allocation ind_buf is failed")
      ind_buf=0

      next_i=0
      label_i: do i=1,n_ind_st
         m=0
         label_j: do j=1,3
            next_i=next_i+1
            neighbour=1
            i1 = tri_index(i, j)
            ind_buf(next_i,neighbour)=i1

            label_k: do k=1,3
               if (k==j) cycle label_k
               i2 = tri_index(i, k)
               xyz_t_buf=(xyz_vertex(i1,:)+xyz_vertex(i2,:))/2.0_r8_kind
               xyz_tt=(radius/sqrt(dot_product(xyz_t_buf,xyz_t_buf)))*xyz_t_buf

               label_l: do l=1,n_xyz_st
                  if(abs(xyz_tt(1)-xyz_vertex(l,1))<small .and. &
                       abs(xyz_tt(2)-xyz_vertex(l,2))<small .and. &
                       abs(xyz_tt(3)-xyz_vertex(l,3))<small) then
                     neighbour=neighbour+1
                     ind_buf(next_i,neighbour)=l
                     m=m+1
                     new_numbers(m)=l
                     cycle label_k
                  endif
               enddo label_l
               n_xyz_st=n_xyz_st+1
               xyz_vertex(n_xyz_st,:)=xyz_tt
               neighbour=neighbour+1
               ind_buf(next_i,neighbour)=n_xyz_st
               m=m+1
               new_numbers(m)=n_xyz_st
            enddo label_k
         enddo label_j
         next_i=next_i+1
         m1=1
         ind_buf(next_i,m1)=new_numbers(1)
         label_n: do n=2,m
            label_n1: do n1=1,n-1
               if(new_numbers(n) == new_numbers(n1)) cycle label_n
            enddo label_n1
            m1=m1+1
            ind_buf(next_i,m1)=new_numbers(n)
            if (m1==3) cycle label_i
         enddo label_n
      enddo label_i

      n_ind_st=next_i
      tri_index = ind_buf

      deallocate (ind_buf, stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:more_triangles: deallocation of ind_buf is failed")
    end subroutine more_triangles
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine put_points_on_shells()
      real(kind=r8_kind) :: r_xyz(3),dr
      integer(kind=i4_kind) :: N_total,status
      integer(kind=i4_kind) :: i_atom,j_atom,i_shell,i_center,n

      N_total=N_spheres*N_centers_on_sphere

      allocate(pot_on_shell(N_total,n_shells),n_points_on_shell(n_shells),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:put_points_on_shells: allocation of pot_on_shell is failed")

      i_sh:do i_shell=1,n_shells

         n_points_on_shell(i_shell)=0
         i_at:do i_atom=1,N_spheres

            i_cen:do i_center=1,N_centers_on_sphere

               r_xyz=xyz_centers(i_center,:)*r_sphere(i_atom,i_shell)/radius+xyz_sphere(i_atom,:)

               j_at:do j_atom=1,N_spheres
                  if(j_atom==i_atom) cycle j_at
                  dr=sqrt(dot_product(r_xyz-xyz_sphere(j_atom,:),r_xyz-xyz_sphere(j_atom,:)))
                  if(dr < r_sphere(j_atom,i_shell)) cycle  i_cen
               end do j_at

               n_points_on_shell(i_shell)=n_points_on_shell(i_shell)+1
               n=n_points_on_shell(i_shell)
               pot_on_shell(n,i_shell)%r=r_xyz
            end do i_cen
         end do i_at
      end do i_sh

      deallocate(xyz_centers,stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module:put_points_on_shells: deallocation of xyz_centers is failed")

    end subroutine put_points_on_shells
    !------------------------------------------------------------

    !------------------------------------------------------------
    subroutine symm_sorted_grid()
      !------------ Modules used --------------------------------------
      use group_module, only : ylm_trafos,sub_group,group_coset, &
           symm_transformation_int,group_num_el,group_coset_decomp
      !------------ Declaration of local variables ------------------
      integer(kind=i4_kind) :: n_equal,n_equal_check,N_total
      real(kind=r8_kind)    :: position(3),position2(3),position3(3)
      type(sub_group) :: local_groups
      type(group_coset) :: cosets
      type(symm_transformation_int) :: point_trafos
      logical, allocatable :: help_dim(:,:)
      type(spacepoint_type),allocatable :: buffer(:,:)
      integer(kind=i4_kind),allocatable :: n_size(:)

      real(kind=r8_kind),parameter :: small = 1.0e-10_r8_kind
      integer (i4_kind) :: i, n, j, k, status, i1, j1, k1, n_p_s
      !------------ Executable code ---------------------------------

      N_total=N_spheres*N_centers_on_sphere
      allocate(help_dim(N_total,n_shells),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: allocation of help_dim is failed")
      help_dim=.false.

      allocate(buffer(N_total,n_shells),n_size(n_shells),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_modulr: allocation of buffer is failed")

      i1_shell:do i1=1,n_shells
         n_p_s=n_points_on_shell(i1)
         n_size(i1)=0
         i=0
         j1_n_p_s:do j1=1,n_p_s
            i=i+1
            help_d: if(.not.help_dim(i,i1)) then
               n_size(i1)=n_size(i1)+1
               ! reorder coordinates of points as
               ! (x,y,z) --> (x,z,y) in order to comply with the
               ! convention for angular momentum l=1
               position(1) = pot_on_shell(j1,i1)%r(1)
               position(2) = pot_on_shell(j1,i1)%r(3)
               position(3) = pot_on_shell(j1,i1)%r(2)
               !
               ! determine local symmetry groups
               !
               ! now apply all symmetry operations to the position of the
               ! surface point
               n = 0
               do j=1,group_num_el
                  position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
                  if (dot_product(position2-position,position2-position) < small) then
                     n = n+1
                  endif
               enddo

               ! allocate group elements
               local_groups%num_el = n
               allocate(local_groups%elements(n))
               if ( status /= 0) call error_handler( &
                    "points_on_cavity_surface: allocation local_groups  is failed")

               ! fill up group elements
               n = 0
               do j=1,group_num_el
                  position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
                  if (dot_product(position2-position,position2-position) < small) then
                     n = n+1
                     local_groups%elements(n) = j
                  end if
               enddo
               !
               ! now determine symmetry equivalent atoms
               !
               call group_coset_decomp(n_equal,local_groups,&
                    cosets,point_trafos%matrix)
               !
               ! search of positions of equal atoms
               !
               allocate(buffer(n_size(i1),i1)%position(3,n_equal), stat=status)
               if (status .ne. 0) call error_handler( &
                    "potential_calc_module: allocation of buffer%position is failed")
               n_equal_check=0
               j_n_equal: do j=1,n_equal
                  position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets%elements(1,j)),position)
                  position3(1) = position2(1)
                  position3(2) = position2(3)
                  position3(3) = position2(2)
                  k=0
                  k1_n_p_s:do k1=1,n_p_s
                     k=k+1
                     if(.not.help_dim(k,i1)) then
                        if (sqrt(dot_product(position3-pot_on_shell(k1,i1)%r, &
                             position3-pot_on_shell(k1,i1)%r)) <= small) then
                           help_dim(k,i1)=.true.
                           n_equal_check=n_equal_check+1
                           buffer(n_size(i1),i1)%position(:,j)=pot_on_shell(k1,i1)%r
                           buffer(n_size(i1),i1)%n_equal_points=n_equal
                        endif
                     end if
                  end do k1_n_p_s
               end do j_n_equal
               if (n_equal_check /= n_equal) call error_handler( &
                    "potential_calc_module: The program has calculated more equivalent points than &
                    & it had been able to find on the shell. Who fib from us? I do not know")

               deallocate(local_groups%elements, &
                    point_trafos%matrix,cosets%elements,stat=status)
               if (status .ne. 0 ) call error_handler( &
                    "potential_calc_module: deallocation of points helpers failed")
            end if help_d
         end do j1_n_p_s
      end do i1_shell

      N_points=sum(n_size)

      allocate(point_in_space(N_points),stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: allocation of point_in_space is failed")

      j=0
      do j1=1,n_shells
         do i=1,n_size(j1)
            j=j+1
            point_in_space(j)%n_equal_points=buffer(i,j1)%n_equal_points

            allocate(point_in_space(j)%position(3,buffer(i,j1)%n_equal_points), stat=status)
            if (status .ne. 0 ) call error_handler( &
                 "potential_calc_module: allocation of point_in_space(j)%position is failed")

            point_in_space(j)%position=buffer(i,j1)%position

            deallocate(buffer(i,j1)%position, stat=status)
            if (status .ne. 0 ) call error_handler( &
                 "potential_calc_module: deallocation of buffer(i,j1)%position is failed")
         end do
      end do

      deallocate(buffer,n_size,n_points_on_shell,stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: deallocation of buffer is failed")
      deallocate(help_dim,stat=status)
      if ( status /= 0) call error_handler( &
           "potential_calc_module: deallocation of help_dim  is failed")

       deallocate(pot_on_shell, stat=status)
       if ( status /= 0) call error_handler( &
              "potential_calc_module: deallocation pot_on_shell is failed")

    end subroutine symm_sorted_grid

  end subroutine calc_shell_grid
!********************************************************************

!********************************************************************
  subroutine collect_poten_3d()
    !
    ! Runs on all workers.
    !
    use density_data_module, only: density_data_free
    !== End of interface ============================================


    call start_read_poten_e()
    call get_poten_n()
!!$    call get_poten_pc() !@@@@@@@@@@@@@@@@@@@@@@@

    call bounds_free_poten()
    call density_data_free()
  end subroutine collect_poten_3d
!********************************************************************

!********************************************************************
  subroutine calc_poten_derive_charges()
    !------------ Modules used --------------------------------------
    use occupation_module, only : get_charge
    use unique_atom_module, only : N_unique_atoms,unique_atoms
!!$    use pointcharge_module, only: pointcharge_N, pointcharge_array !@@@@@@@@@@@@@@
    use dipole_module, only : dipole_total
    !------------ Declaration of formal parameters ------------------
    !== End of interface ============================================
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: N_charges,N_dim,N_dim_old
    integer(kind=i4_kind), allocatable :: IPIV(:),b_help(:)
    real(kind=r8_kind), allocatable :: V_fix_at(:),a_matrix(:,:),b_vector(:,:)
    real(kind=r8_kind), allocatable :: b_buffer(:,:), a_buffer(:,:)
    real(kind=r8_kind) :: dipmom_fix_at(3),rp(3),ra0(3),ra1(3),dr0,dr1
    integer(kind=i4_kind) :: n_eq_at0,n_eq_at1,n_eq_p,n_constr,n_symm,n_symms(30),s_max
    integer(kind=i4_kind) :: l_max(1), numbers(100,30),N_cycles,N_fixed,total_n(200)
    logical :: use_constr(3)
    real(kind=r8_kind),parameter :: small = 1.0e-10_r8_kind
    real(kind=r8_kind),parameter :: a_strong=0.001_r8_kind
    real(kind=r8_kind),parameter :: a_weak=0.0005_r8_kind
    integer(kind=i4_kind) :: i,j,k,l,m,n,ii,jj,i1,j1,status,info,icycle
    integer (i4_kind) :: is, ks, ls
!!$    real(r8_kind) :: pc_charge !@@@@@@@@@@@@@
    !------------ Executable code -----------------------------------

    N_cycles=1
    if(symm_restr .and. n_fixed_atoms > 0) N_cycles=2
    N_fixed=0

    N_cyc: do icycle=1,N_cycles

       if(symm_restr .and. icycle == 1) then
          N_fixed=n_fixed_atoms
          n_fixed_atoms=0
       end if

       allocate(V_fix_at(N_points),stat=status)
       if ( status /= 0) call error_handler( &
            "potential_calc_module: allocation of V_fix_at is failed")
       V_fix_at=0.0_r8_kind

       if(n_fixed_atoms > 0) then
          do i=1,N_points
             rp=point_in_space(i)%position(:,1)
             do j=1,n_fixed_atoms
                k=fixed_atoms(j)
                if(k <= N_unique_atoms) n_eq_at0=unique_atoms(k)%N_equal_atoms
!!$                if(k > N_unique_atoms) &                                          !@@@@@@@@
!!$                     n_eq_at0=pointcharge_array(k-N_unique_atoms)%N_equal_charges !@@@@@@@@
                do l=1,n_eq_at0
                   if(k <= N_unique_atoms) ra0=unique_atoms(k)%position(:,l)
!!$                   if(k > N_unique_atoms) &                                       !@@@@@@@@
!!$                        ra0=pointcharge_array(k-N_unique_atoms)%position(:,l)     !@@@@@@@@
                   dr0=sqrt(dot_product(ra0-rp,ra0-rp))
                   V_fix_at(i)=V_fix_at(i)+fixed_charges(j)/dr0
                end do
             end do
          end do

          dipmom_fix_at=0.0_r8_kind
          if(n_fixed_atoms > 0) then
             do j=1,n_fixed_atoms
                k=fixed_atoms(j)
                if(k <= N_unique_atoms) n_eq_at0=unique_atoms(k)%N_equal_atoms
!!$                if(k > N_unique_atoms) &                                          !@@@@@@@@
!!$                     n_eq_at0=pointcharge_array(k-N_unique_atoms)%N_equal_charges !@@@@@@@@
                do l=1,n_eq_at0
                   if(k <= N_unique_atoms) ra0=unique_atoms(k)%position(:,l)
!!$                   if(k > N_unique_atoms) &                                       !@@@@@@@@
!!$                        ra0=pointcharge_array(k-N_unique_atoms)%position(:,l)     !@@@@@@@@
                   dipmom_fix_at=dipmom_fix_at+fixed_charges(j)*ra0
                end do
             end do
          end if
       end if

       N_charges = 0
       lab1:do i=1,N_unique_atoms
          if(n_fixed_atoms > 0) then
             do j=1,n_fixed_atoms
                if(i==fixed_atoms(j)) cycle lab1
             end do
          end if
          N_charges = N_charges + unique_atoms(i)%N_equal_atoms
       end do lab1
!!$       labp1: do i=1,pointcharge_N                                       !@@@@@@@@@@@@@@@@
!!$          if(n_fixed_atoms > 0) then                                     !@@@@@@@@@@@@@@@@
!!$             do j=1,n_fixed_atoms                                        !@@@@@@@@@@@@@@@@
!!$                if(i+N_unique_atoms==fixed_atoms(j)) cycle labp1         !@@@@@@@@@@@@@@@@
!!$             end do                                                      !@@@@@@@@@@@@@@@@
!!$          end if                                                         !@@@@@@@@@@@@@@@@
!!$          N_charges = N_charges + pointcharge_array(i)%N_equal_charges   !@@@@@@@@@@@@@@@@
!!$       end do labp1                                                      !@@@@@@@@@@@@@@@@

       if(charge_constr) n_constr=1
       if(dipole_constr) then
          use_constr=.true.
          n_constr=4
          lab2:do i=1,N_unique_atoms
             if(n_fixed_atoms > 0) then
                do j=1,n_fixed_atoms
                   if(i==fixed_atoms(j)) cycle lab2
                end do
             end if
             n_eq_at0=unique_atoms(i)%N_equal_atoms
             do j=1,n_eq_at0
                ra0=unique_atoms(i)%position(:,j)
                do k=1,3
                   if(abs(ra0(k)) <= small) then
                      if(use_constr(k)) n_constr=n_constr-1
                      use_constr(k)=.false.
                   end if
                end do
             end do
          end do lab2
!!$          labp2:do i=1,pointcharge_N                                !@@@@@@@@@@@@@@@
!!$             if(n_fixed_atoms > 0) then                             !@@@@@@@@@@@@@@@
!!$                do j=1,n_fixed_atoms                                !@@@@@@@@@
!!$                   if(i+N_unique_atoms==fixed_atoms(j)) cycle labp2 !@@@@@@@@@@@@
!!$                end do                                              !@@@@@@@@@@@@
!!$             end if                                                 !@@@@@@@@@@@@@
!!$             n_eq_at0=pointcharge_array(i)%N_equal_charges          !@@@@@@@@@@@@@@@
!!$             do j=1,n_eq_at0                                        !@@@@@@@@@@@@@@@
!!$                ra0=pointcharge_array(i)%position(:,j)              !@@@@@@@@@@@@@@@
!!$                do k=1,3                                            !@@@@@@@@@@@@@@@@
!!$                   if(abs(ra0(k)) <= small) then                    !@@@@@@@@@@@@@@@
!!$                      if(use_constr(k)) n_constr=n_constr-1         !@@@@@@@@@@@@@@@
!!$                      use_constr(k)=.false.                         !@@@@@@@@@@@@@@@
!!$                   end if                                           !@@@@@@@@@@@@@@@
!!$                end do                                              !@@@@@@@@@@@@@@@@
!!$             end do                                                 !@@@@@@@@@@@@@@@
!!$          end do labp2                                              !@@@@@@@@@@@@@@@@
       end if
       N_dim=N_charges+n_constr
       if(n_constr > N_charges) then
          dipole_constr = .false.
          charge_constr = .true.
          N_dim=N_charges+1
       end if
       N_dim_old=N_dim

       allocate(a_matrix(N_dim,N_dim),b_vector(N_dim,1),b_help(N_dim),stat=status)
       if ( status /= 0) call error_handler( &
            "potential_calc_module: allocation of a_matrix or b_vector is failed")
       a_matrix=0.0_r8_kind
       b_vector=0.0_r8_kind

       ii=0
       lab_i:do i=1,N_unique_atoms
!!$       lab_i:do i=1,N_unique_atoms+pointcharge_N        !@@@@@@@@@@@@@@
          if(n_fixed_atoms > 0) then
             do j=1,n_fixed_atoms
                if(i==fixed_atoms(j)) cycle lab_i
             end do
          end if

          if(i <= N_unique_atoms) n_eq_at0=unique_atoms(i)%N_equal_atoms
!!$          if(i > N_unique_atoms) n_eq_at0=pointcharge_array(i-N_unique_atoms)%N_equal_charges !@@@@@@@@@@@@@@
          lab_j:do j=1,n_eq_at0

             ii=ii+1
             if(i <= N_unique_atoms) ra0=unique_atoms(i)%position(:,j)
!!$             if(i > N_unique_atoms) ra0=pointcharge_array(i-N_unique_atoms)%position(:,j) !@@@@@@@@@@@@@
             b_vector(ii,1)=0.0_r8_kind
             lab_i1:do i1=1,N_points

                n_eq_p=point_in_space(i1)%n_equal_points
                lab_j1:do j1=1,n_eq_p

                   rp=point_in_space(i1)%position(:,j1)
                   dr0=sqrt(dot_product(ra0-rp,ra0-rp))
                   b_vector(ii,1)=b_vector(ii,1)+ &
                        (V_pot_e(i1)+V_pot_n(i1)-n_eq_p*V_fix_at(i1))/(n_eq_p*dr0)
!!$                        (V_pot_e(i1)+V_pot_n(i1)+V_pot_pc(i1)-n_eq_p*V_fix_at(i1))/(n_eq_p*dr0)  !@@@@@@@@@@@
                end do lab_j1

             end do lab_i1

             jj=0
             lab_k:do k=1,N_unique_atoms
!!$             lab_k:do k=1,N_unique_atoms+pointcharge_N        !@@@@@@@@@@@@@@
                if(n_fixed_atoms > 0) then
                   do l=1,n_fixed_atoms
                      if(k==fixed_atoms(l)) cycle lab_k
                   end do
                end if

                if(k <= N_unique_atoms) n_eq_at1=unique_atoms(k)%N_equal_atoms
!!$                if(k > N_unique_atoms) n_eq_at1=pointcharge_array(k-N_unique_atoms)%N_equal_charges !@@@@@@@@@@@@@@
                lab_l:do l=1,n_eq_at1

                   jj=jj+1
                   if(k <= N_unique_atoms) ra1=unique_atoms(k)%position(:,l)
!!$                   if(k > N_unique_atoms) ra1=pointcharge_array(k-N_unique_atoms)%position(:,l) !@@@@@@@@@@@@@
                   labb_i1:do i1=1,N_points

                      n_eq_p=point_in_space(i1)%n_equal_points
                      labb_j1:do j1=1,n_eq_p

                         rp=point_in_space(i1)%position(:,j1)
                         dr0=sqrt(dot_product(ra0-rp,ra0-rp))
                         dr1=sqrt(dot_product(ra1-rp,ra1-rp))
                         a_matrix(ii,jj)=a_matrix(ii,jj)+1.0_r8_kind/(dr0*dr1) !+2.0_r8_kind !???AS
                      end do labb_j1

                   end do labb_i1

                end do lab_l

             end do lab_k

          end do lab_j

       end do lab_i

       if(charge_constr) m=1
       if(dipole_constr) m=4
       ii=0
       lab3:do i=1,N_unique_atoms
!!$       lab3:do i=1,N_unique_atoms+pointcharge_N       !@@@@@@@@@@@@@@
          if(n_fixed_atoms > 0) then
             do j=1,n_fixed_atoms
                if(i==fixed_atoms(j)) cycle lab3
             end do
          end if

          if(i <= N_unique_atoms) n_eq_at0=unique_atoms(i)%N_equal_atoms
!!$          if(i > N_unique_atoms) n_eq_at0=pointcharge_array(i-N_unique_atoms)%N_equal_charges !@@@@@@@@@@@@@@
          do k=1,n_eq_at0
             ii=ii+1
             if(i <= N_unique_atoms) ra0=unique_atoms(i)%position(:,k)
!!$             if(i > N_unique_atoms) ra0=pointcharge_array(i-N_unique_atoms)%position(:,k) !@@@@@@@@@@@@@
             l=0
             do j=1,m
                if(j == 1) then
                   l=l+1
                   a_matrix(ii,N_charges+l)=1.0_r8_kind
                end if
                if(j > 1) then
                   if(.not.use_constr(j-1)) cycle
                   l=l+1
                   a_matrix(ii,N_charges+l)=ra0(j-1)
                end if
                a_matrix(N_charges+l,ii)= a_matrix(ii,N_charges+l)
             end do
          end do
       end do lab3

!!$       pc_charge=0.0_r8_kind                            !@@@@@@@@@@@@
!!$       do i=1,pointcharge_N                             !@@@@@@@@@@@@
!!$          do j=1,pointcharge_array(i)%N_equal_charges   !@@@@@@@@@@@
!!$             pc_charge=pc_charge+pointcharge_array(i)%Z !@@@@@@@@@@@
!!$          end do                                        !@@@@@@@@@@@@@
!!$       end do                                           !@@@@@@@@@@@@@2

       l=0
       do i=1,m
          if(i == 1) then
             l=l+1
             b_vector(N_charges+l,1)=get_charge()
!!$             b_vector(N_charges+l,1)=get_charge()+pc_charge !@@@@@@@@@@@@@
             if(n_fixed_atoms > 0) then
                do j=1,n_fixed_atoms
                   k=fixed_atoms(j)
                   if(k <= N_unique_atoms) n_eq_at0=unique_atoms(k)%N_equal_atoms
!!$                   if(k > N_unique_atoms) n_eq_at0=pointcharge_array(k-N_unique_atoms)%N_equal_charges !@@@@@@@@@@@@@@
                   b_vector(N_charges+l,1)=b_vector(N_charges+l,1)-n_eq_at0*fixed_charges(j)
                end do
             end if
          end if
          if(i > 1) then
             if(.not.use_constr(i-1)) cycle
             l=l+1
             b_vector(N_charges+l,1)=dipole_total(i-1)-dipmom_fix_at(i-1)
          end if
       end do

       do i=1,N_dim_old
          b_help(i)=i
       end do

       if((symm_restr .and. N_cycles==1) .or. (symm_restr .and. icycle==2)) then
          !resorting symm_atoms list to garantee ascending order of
          !unique atomic numbers
          ks=0
          do is=1,n_symm_groups
             do i=1,n_symm_atoms(is)-1
                k=n_symm_atoms(is)-i+1
                l_max=maxloc(symm_atoms(ks+1:ks+k))
                s_max=symm_atoms(ks+l_max(1))
                symm_atoms(l_max(1)+ks:ks+k)=eoshift(symm_atoms(l_max(1)+ks:ks+k),shift=1)
                symm_atoms(ks+k)=s_max
             end do
             ks=ks+n_symm_atoms(is)
          end do

          !Due to forced equivalence of some atoms sizes of a_matrix
          !and b_vector have to be reduced by summing up equivalent
          !rows and columns
          l=0; ks=0
          do is=1,n_symm_groups
             ls=0
             do i=1,n_symm_atoms(is)
                j=symm_atoms(i+ks)
                m=0
                aaa:do k=1,j-1
                   do i1=1,n_fixed_atoms
                      if(fixed_atoms(i1)==k) cycle aaa
                   end do
                   if(k <= N_unique_atoms) m=m+unique_atoms(k)%N_equal_atoms
!!$                   if(k > N_unique_atoms) m=m+pointcharge_array(k-N_unique_atoms)%N_equal_charges !@@@@@@@@@@@@@@
                end do aaa
                if(j <= N_unique_atoms) n=unique_atoms(j)%N_equal_atoms
!!$                if(j > N_unique_atoms) n=n+pointcharge_array(j-N_unique_atoms)%N_equal_charges !@@@@@@@@@@@@@@
                do k=m+1,m+n
                   ls=ls+1
                   l=l+1
                   numbers(ls,is)=k
                end do
             end do
             ks=ks+n_symm_atoms(is)
             n_symms(is)=ls
          end do

          n_symm=l
          !summing up columns
          ks=0
          do is=1,n_symm_groups
             i=numbers(1,is)
             do j=2,n_symms(is)
                k=numbers(j,is)
                a_matrix(:,i)=a_matrix(:,i)+a_matrix(:,k)
                b_help(k)=i
             end do
             !summing up rows
             do j=2,n_symms(is)
                k=numbers(j,is)
                a_matrix(i,:)=a_matrix(i,:)+a_matrix(k,:)
                b_vector(i,1)=b_vector(i,1)+b_vector(k,1)
             end do
             total_n(ks+1:ks+n_symms(is)-1)=numbers(2:n_symms(is),is)
             ks=n_symms(is)-1
          end do

          ks=n_symm-n_symm_groups
          do i=1,ks-1
             l_max=maxloc(total_n(1:ks-i+1))
             s_max=total_n(l_max(1))
             total_n(l_max(1):ks-i+1)=eoshift(total_n(l_max(1):ks-i+1),shift=1)
             total_n(ks-i+1)=s_max
          end do

          do i=ks,1,-1
             j=total_n(i)
             !erasing columns
             a_matrix(:,j:N_dim)=eoshift(a_matrix(:,j:N_dim),shift=1,dim=2)
             do l=1,N_dim_old-1
                if(b_help(l) > j) b_help(l)=b_help(l)-1
             end do
          end do
          do i=ks,1,-1
             j=total_n(i)
             !erasing rows
             a_matrix(j:N_dim,:)=eoshift(a_matrix(j:N_dim,:),shift=1,dim=1)
             b_vector(j:N_dim,1)=eoshift(b_vector(j:N_dim,1),shift=1,dim=1)
          end do

          N_dim=N_dim-n_symm+n_symm_groups
       end if

       allocate(b_buffer(N_dim,1),a_buffer(N_dim,N_dim),stat=status)
       if ( status /= 0) call error_handler( &
            "potential_calc_module: allocation b_buffer is failed")
       b_buffer(:,1)=b_vector(1:N_dim,1)
       a_buffer=a_matrix(1:N_dim,1:N_dim)

       allocate(IPIV(N_dim),stat=status)
       if(status /= 0) call error_handler("potential_calc_module: allocation ipiv is failed")

       call dgesv(N_dim,1,a_buffer,N_dim,IPIV,b_buffer,N_dim,info)
       !now b_bufferq contains charges

       if(resp_restr) then
          if(symm_restr .and. icycle==2) resp_a=a_strong
          call nonlinear_solver()
       end if

       b_vector(1:N_dim,1)=b_buffer(:,1)

       if(.not.symm_restr .or. icycle == 2) call print_results()
       if(symm_restr .and. N_cycles == 1) call print_results()

       if(symm_restr .and. N_fixed > 0) then
          n_fixed_atoms=N_fixed
          N_fixed=0
          do i=1,n_fixed_atoms
             j=fixed_atoms(i)
             l=0
             do k=1,j-1
                if(k <= N_unique_atoms) l=l+unique_atoms(k)%N_equal_atoms
!!$                if(k > N_unique_atoms) l=l+pointcharge_array(k-N_unique_atoms)%N_equal_charges !@@@@@@@@@@@@@@
             end do
             l=l+1
             fixed_charges(i)=b_vector(l,1)
          end do
       end if

       deallocate(V_fix_at,stat=status)
       if ( status /= 0) call error_handler( &
            "potential_calc_module: deallocation of v_fix_at is failed")

       deallocate(a_matrix,b_vector,b_help,IPIV,stat=status)
       if ( status /= 0) call error_handler( &
            "potential_calc_module: deallocation of a_matrix or b_vector is failed")

       deallocate(b_buffer,a_buffer,stat=status)
       if ( status /= 0) call error_handler( &
            "potential_calc_module: deallocation b_buffer is failed")

       if((.not.symm_restr .or. icycle == 2) .or. &
            (symm_restr .and. N_cycles == 1)) then
          call dealloc_space_points()
          call deallocate_pot()
       end if

    end do N_cyc

  contains
    !---------------------------------------------------------
    subroutine nonlinear_solver

        integer (i4_kind) :: in, jn, statusn
        integer(i4_kind), parameter :: max_iter=100
        real(r8_kind), parameter :: qtol=0.1e-6_r8_kind
        real(r8_kind), allocatable :: qr(:)
        real(r8_kind) :: qdiff

        allocate(qr(N_dim),stat=statusn)
        if(statusn /= 0) call error_handler("potential_calc_module: allocation qr is failed")

        qr=b_buffer(:,1)

        main: do in=1,max_iter
           a_buffer=a_matrix(1:N_dim,1:N_dim)
           b_buffer(:,1)=b_vector(1:N_dim,1)

           !do restraint correction
           do jn=1,N_dim-1
              a_buffer(jn,jn)=a_buffer(jn,jn)+resp_a/sqrt(qr(jn)*qr(jn)+resp_b*resp_b)
              if(a_buffer(jn,jn) <= 1.0e-10_r8_kind) a_buffer(jn,jn) = 1.0e-10_r8_kind
           end do

           call dgesv(N_dim,1,a_buffer,N_dim,IPIV,b_buffer,N_dim,info)
           !now b_vector contains charges

           !check convergence
           qdiff=sqrt(dot_product(qr(1:N_dim-1)-b_buffer(1:N_dim-1,1), &
                qr(1:N_dim-1)-b_buffer(1:N_dim-1,1)))
           if(qdiff <= qtol) exit main

           qr=b_buffer(:,1)

        end do main

        deallocate(qr,stat=statusn)
        if(statusn /= 0) call error_handler("potential_calc_module: deallocation qr is failed")

      end subroutine nonlinear_solver
      !---------------------------------------------------------

      !---------------------------------------------------------
      real(kind=r8_kind) function calc_esp_pdc(rp)
        !------------ Declaration of local variables --------------------
        real(kind=r8_kind) :: ra(3),rp(3),dr,Q
        integer (i4_kind) :: i, j, ii, jj, kk
        !------------ Executable code -----------------------------------
        calc_esp_pdc=0.0_r8_kind
        ii=0; kk=0
        lab1:do i=1,N_unique_atoms
           if(n_fixed_atoms > 0) then
              do j=1,n_fixed_atoms
                 if(i==fixed_atoms(j)) cycle lab1
              end do
           end if
           do j=1,unique_atoms(i)%N_equal_atoms
              ii=ii+1
              jj=b_help(ii)
              Q=b_vector(jj,1)
              ra=unique_atoms(i)%position(:,j)
              dr=sqrt(dot_product(ra-rp,ra-rp))
              calc_esp_pdc=calc_esp_pdc+Q/dr
           end do
        end do lab1
!!$        labp1:do i=1,pointcharge_N                                 !@@@@@@@@@@@@@@
!!$           if(n_fixed_atoms > 0) then                              !@@@@@@@@@@@@@
!!$              do j=1,n_fixed_atoms                                 !@@@@@@@@@@@@@@
!!$                 if(i+N_unique_atoms==fixed_atoms(j)) cycle labp1  !@@@@@@@@@@@@@
!!$              end do                                               !@@@@@@@@@@@@@
!!$           end if                                                  !@@@@@@@@@@@@@@@
!!$           do j=1,pointcharge_array(i)%N_equal_charges             !@@@@@@@@@@@@@@@
!!$              ii=ii+1                                              !@@@@@@@@@@@@@@@
!!$              jj=b_help(ii)                                        !@@@@@@@@@@@@@@@
!!$              Q=b_vector(jj,1)                                     !@@@@@@@@@@@@@@@
!!$              ra=pointcharge_array(i)%position(:,j)                !@@@@@@@@@@@@@@@@
!!$              dr=sqrt(dot_product(ra-rp,ra-rp))                    !@@@@@@@@@@@@@@@@
!!$              calc_esp_pdc=calc_esp_pdc+Q/dr                       !@@@@@@@@@@@@@@@@
!!$           end do                                                  !@@@@@@@@@@@@@@@@@
!!$        end do labp1                                               !@@@@@@@@@@@@@@@@@
      end function calc_esp_pdc
      !---------------------------------------------------------

      !---------------------------------------------------------
      subroutine print_results()
        !------------ Modules used --------------------------------------
#ifdef WITH_MOLMECH
        use operations_module, only : operations_qm_mm_new
#endif
        !------------ Declaration of local variables --------------------
        real(kind=r8_kind) :: esp_qm,desp_qm,esp_pdc,desp,dipmom_pdc(3)
        integer (i4_kind) :: i, j, k, ii, nn, jj, kk
        logical :: save_pdc
        integer :: pdc_unit
        !------------ Executable code -----------------------------------

        save_pdc=.false.
#ifdef WITH_MOLMECH
        if(operations_qm_mm_new) then
           save_pdc=.true.
           pdc_unit=openget_iounit(trim(inpfile('pdc.save')),  &
                form='formatted', status='unknown')
        end if
#endif
        write(output_unit,*) '================= Potential Derived Charges =================='
        write(output_unit,*) '------------------Merz-Kollman method (ESP) ------------------'
        write(output_unit,*) 'B.H.Besler, K.M.Merz, P.A.Kollman, J.Comp.Chem., 11 (1990) 431'
        write(output_unit,*) 'C.I.Bayly, P.Cieplak, W.D.Cornell, P.A.Kollman,  J.Phys.Chem.,'
        write(output_unit,*) '97 (1993) 10269'
        write(output_unit,*) '--------------------------------------------------------------'
        write(output_unit,*) ''
        if(pdc_out_type >= 2) &
             write(output_unit,'(3x,i4,a48)') &
             n_shells,' shells has been built around molecule (cluster)'
        if(pdc_out_type >= 2) write(output_unit,*) 'to calculate electrostatic potential'
        if(pdc_out_type >= 2) write(output_unit,*) ''

        if(pdc_out_type >= 2) write(output_unit,*) ' Radii of the atomic spheres'
        if(pdc_out_type >= 2) write(output_unit,*) '-----------------------------'
        if(pdc_out_type >= 2) write(output_unit,*) ' Atom        Radius(angstrom)'
        if(pdc_out_type >= 2) write(output_unit,*) '-----------------------------'
        k=1
        do i=1,N_unique_atoms
           if(pdc_out_type >= 2) write(output_unit,'(3x,a4,10x,f10.5)') &
                unique_atoms(i)%name,r_sphere(k,1)*au2ang
           do j=2,n_shells
              if(pdc_out_type >= 2) write(output_unit,'(17x,f10.5)') r_sphere(k,j)*au2ang
           end do
           k=k+unique_atoms(i)%N_equal_atoms
        end do
!!$        do i=1,pointcharge_N                                                                !@@@@@@@@@@
!!$           if(pdc_out_type >= 2) write(output_unit,'(3x,a4,10x,f10.5)') &                   !@@@@@@@@@@@
!!$                pointcharge_array(i)%name,r_sphere(k,1)*au2ang                              !@@@@@@@@@@@
!!$           do j=2,n_shells                                                                  !@@@@@@@@@@@
!!$              if(pdc_out_type >= 2) write(output_unit,'(17x,f10.5)') r_sphere(k,j)*au2ang   !@@@@@@@@@@@@@
!!$           end do                                                                           !@@@@@@@@@@@@@
!!$           k=k+pointcharge_array(i)%N_equal_charges                                         !@@@@@@@@@@@@
!!$        end do                                                                              !@@@@@@@@@@@@2
        deallocate(r_sphere)
        if(pdc_out_type >= 2) write(output_unit,*) '-----------------------------'

        if(pdc_out_type >= 2 .and. poligon_type == 1)  &
             write (output_unit,*) 'The DODECAHEDRON has been inscribed into each sphere'
        if(pdc_out_type >= 2 .and. poligon_type == 2) &
             write (output_unit,*) 'The CUBE has been inscribed into each sphere'
        if(pdc_out_type >= 2 .and. poligon_type == 3) &
             write (output_unit,*) 'The DOBLEPIRAMIDE has been inscribed into each sphere'
        if(pdc_out_type >= 2) &
             write(output_unit,'(a16,i4,a24)') &
             ' On each sphere ',N_centers_on_sphere,' points has been imposed'
        ii=0
        do i=1,N_points
           ii=ii+point_in_space(i)%n_equal_points
        end do
        if(pdc_out_type >= 2) write(output_unit,'(a39,i5)') ' Total number of points on the shells: ',ii
        if(pdc_out_type >= 2) write(output_unit,'(a31,i5)') ' Number of symmetrized points: ',N_points
        if(pdc_out_type >= 2) write(output_unit,*) ''

        if(pdc_out_type == 3) write(output_unit,*) '      Point on shell*           ESP(QM)*   ESP(PDC)*'
        if(pdc_out_type == 3) write(output_unit,*) '    x        y        z '
        if(pdc_out_type == 3) write(output_unit,*) '----------------------------------------------------'
        desp=0.0_r8_kind
        desp_qm=0.0_r8_kind
        ii=0
        do i=1,N_points
           do j=1,point_in_space(i)%n_equal_points
              nn=point_in_space(i)%n_equal_points
              esp_qm=(V_pot_e(i)+V_pot_n(i))/nn !-V_fix_at(i)
!!$              esp_qm=(V_pot_e(i)+V_pot_n(i)+V_pot_pc(i))/nn       !@@@@@@@@@@@@@@@
              desp_qm=desp_qm+(esp_qm)**2
              esp_pdc=calc_esp_pdc(point_in_space(i)%position(:,j))+V_fix_at(i)
              desp=desp+(esp_qm-esp_pdc)**2
              if(pdc_out_type == 3) &
                   write(output_unit,'(3(1x,f8.4),2(4x,f8.4))') &
                   point_in_space(i)%position(:,j),esp_qm,esp_pdc
           end do
        end do
        if(pdc_out_type == 3) write(output_unit,*) '----------------------------------------------------'
        if(pdc_out_type == 3) write(output_unit,*) ' *) all data in au'
        if(pdc_out_type == 3) write(output_unit,*) ''

        if(pdc_out_type >= 1) write(output_unit,*) ' Unique number    Atom       Charge'
        if(pdc_out_type >= 1) write(output_unit,*) '-------------------------------------'

        ii=0; kk=0
        lab1:do i=1,N_unique_atoms
!!$        lab1:do i=1,N_unique_atoms+pointcharge_N          !@@@@@@@@@@@@@@@@@@@
           if(n_fixed_atoms > 0) then
              do j=1,n_fixed_atoms
                 if(i==fixed_atoms(j)) then
                    if(i <= N_unique_atoms) then
                       do k=1,unique_atoms(i)%N_equal_atoms
                          if(pdc_out_type >= 1)write(output_unit,'(6x,i3,11x,a4,3x,f10.5,a11)') &
                               i,unique_atoms(i)%name,fixed_charges(j),' predefined'
                       end do
                    end if
!!$                    if(i > N_unique_atoms) then                                                          !@@@@@@@@@
!!$                       do k=1,pointcharge_array(i-N_unique_atoms)%N_equal_charges                        !@@@@@@@@
!!$                          if(pdc_out_type >= 1)write(output_unit,'(6x,i3,11x,a4,3x,f10.5,a11)') &        !@@@@@@@@@@
!!$                               i,pointcharge_array(i-N_unique_atoms)%name,fixed_charges(j),' predefined' !@@@@@@@@@@
!!$                       end do                                                                            !@@@@@@@@@
!!$                    end if                                                                               !@@@@@@@@@@@
                    cycle lab1
                 end if
              end do
           end if
           if(i <= N_unique_atoms) then
              do j=1,unique_atoms(i)%N_equal_atoms
                 ii=ii+1
                 jj=b_help(ii)
                 if(pdc_out_type >= 1) then
                    write (output_unit,'(6x,i3,11x,a4,3x,f10.5)') &
                         i,unique_atoms(i)%name,b_vector(jj,1)
                    if(save_pdc) write (pdc_unit,*) i,b_vector(jj,1)
                 end if
              end do
           end if
!!$           if(i > N_unique_atoms) then                                                                                   !@@@@@
!!$              do j=1,pointcharge_array(i-N_unique_atoms)%N_equal_charges                                                 !@@@@@
!!$                 ii=ii+1                                                                                                 !@@@@@
!!$                 jj=b_help(ii)                                                                                           !@@@@@
!!$                 if(pdc_out_type >= 1) write &                                                                           !@@@@@
!!$                      (output_unit,'(6x,i3,11x,a4,3x,f10.5)') i,pointcharge_array(i-N_unique_atoms)%name,b_vector(jj,1) !@@@@@
!!$              end do                                                                                                     !@@@@@
!!$           end if                                                                                                        !@@@@@
        end do lab1
        if(pdc_out_type >= 1)write(output_unit,*) '-------------------------------------'
        if(pdc_out_type >= 1)write(output_unit,'(a11,f8.5)') ' ESP rrms: ',sqrt(desp/desp_qm)

        ii=0; kk=0
        dipmom_pdc=0.0_r8_kind
        lab2:do i=1,N_unique_atoms
!!$        lab2:do i=1,N_unique_atoms+pointcharge_N          !@@@@@@@@@@@@@@@@@@@
           if(n_fixed_atoms > 0) then
              do j=1,n_fixed_atoms
                 if(i==fixed_atoms(j)) then
                    if(i <= N_unique_atoms) then
                       do k=1,unique_atoms(i)%N_equal_atoms
                          dipmom_pdc=dipmom_pdc+fixed_charges(j)*unique_atoms(i)%position(:,k)
                       end do
                    end if
!!$                    if(i > N_unique_atoms) then                                                                   !@@@@@@@@@
!!$                       do k=1,pointcharge_array(i-N_unique_atoms)%N_equal_charges                                 !@@@@@@@@@
!!$                          dipmom_pdc=dipmom_pdc+fixed_charges(j)*pointcharge_array(i-N_unique_atoms)%position(:,k)!@@@@@@@@@
!!$                       end do                                                                                     !@@@@@@@@@
!!$                    end if                                                                                        !@@@@@@@@@
                    cycle lab2
                 end if
              end do
           end if
           if(i <= N_unique_atoms) then
              do j=1,unique_atoms(i)%N_equal_atoms
                 ii=ii+1
                 jj=b_help(ii)
                 dipmom_pdc=dipmom_pdc+b_vector(jj,1)*unique_atoms(i)%position(:,j)
              end do
           end if
!!$           if(i > N_unique_atoms) then                                                                    !@@@@@@@@@
!!$              do j=1,pointcharge_array(i-N_unique_atoms)%N_equal_charges                                  !@@@@@@@@@
!!$                 ii=ii+1                                                                                  !@@@@@@@@@
!!$                 jj=b_help(ii)                                                                            !@@@@@@@@@
!!$                 dipmom_pdc=dipmom_pdc+b_vector(jj,1)*pointcharge_array(i-N_unique_atoms)%position(:,j)   !@@@@@@@@@
!!$              end do                                                                                      !@@@@@@@@@
!!$           end if                                                                                         !@@@@@@@@@
        end do lab2
        if(pdc_out_type >= 1) write(output_unit,'(a23,3(1x,f8.4))') ' Dipole moment of PDC: ', &
             dipmom_pdc
        if(pdc_out_type >= 1 .and. dipole_constr) &
             write(output_unit,'(a23,3(1x,f8.4))') ' QM dipole moment:     ', &
             dipole_total-dipmom_fix_at

        write(output_unit,*) ''
        write(output_unit,*) '================= Potential Derived Charges =================='

#ifdef WITH_MOLMECH
        if(save_pdc) call returnclose_iounit(pdc_unit)
#endif
      end subroutine print_results

  end subroutine calc_poten_derive_charges

end module potential_calc_module

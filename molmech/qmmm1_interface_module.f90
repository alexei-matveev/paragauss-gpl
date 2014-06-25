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
module qmmm1_interface_module
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module
#ifdef _COMPAC_FORTRAN
  use datatype
#endif
  use iounitadmin_module
  use filename_module
  use common_data_module

  implicit none
  private
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  type, public :: gx_file
     real(r8_kind) :: number
     real(r8_kind) :: x,y,z
     real(r8_kind) :: charge
     integer :: i_uq
     character(len=5) :: name
  end type gx_file
  type(gx_file), allocatable, public :: gx_qmmm(:)

  type, public :: gx_grad
     real(r8_kind) :: x,y,z
  end type gx_grad
  type(gx_grad), allocatable, public :: gx_qmmm_grad(:)

  real(r8_kind), public :: Energy_qmmm, Energy_qm, Energy_mm

  integer, public :: N_gx_atoms, N_gx_qm_atoms
  !------------ public functions and subroutines ------------------
  public read_gx_qmmm,qmmm2pc,QMfield_at_mm_points,qm_grads_to_qmmm1, &
       mm_grads_to_qmmm1,write_gx_qmmm,read_qmmm1,def_qm_mm_1_tasks
  !================================================================
  ! End of public interface of module
  !================================================================
  !------------ Declaration of private constants and variables ----
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine def_qm_mm_1_tasks()
    !------------ modules used ------------------------------------
    use qmmm_interface_module, only : qm_mm_1_task
    use operations_module, only : operations_potential,operations_post_scf, &
         operations_gradients,operations_geo_opt,operations_solvation_effect
    use potential_calc_module, only : pdc
    ! --- Declaration of local variables --------------------------
    ! -------------------------------------------------------------

    if(qm_mm_1_task == 0) then
       operations_potential=.false.
    elseif(qm_mm_1_task == 1) then
       operations_potential=.true.
       operations_post_scf =.false.
       operations_gradients=.false.
       operations_geo_opt=.false.
       operations_solvation_effect=.false.
       pdc=.true.
    end if

  end subroutine def_qm_mm_1_tasks
  !****************************************************************

  !****************************************************************
  subroutine read_gx_qmmm()
    !------------ modules used ----------------------------------
    use unique_atom_module
    ! --- Declaration of local variables --------------------
    real(r8_kind) :: iwork
    integer :: gxf,ios,status
    character(len=200) :: buffer
    real(r8_kind) :: coor(4)
    integer :: uniq
    integer :: n_at,i,i_ua,counter_equal

    type read_gx
       type(gx_file) :: data
       type(read_gx), pointer :: next_data
    end type read_gx
    type(read_gx), target :: first_data
    type(read_gx), pointer :: current_data, tmp_data, del
    !------------ Executable code --------------------------------

    gxf = openget_iounit(file=trim(inpfile('gxfile')), status='old', &
         form='formatted')

    current_data=>first_data
    nullify(current_data%next_data)

    n_at=0
    do
       read(gxf,'(a200)',err=10) buffer
       coor=zero
       read(buffer,*,iostat=ios) coor(1:4),uniq
       if(ios /= 0) then
          coor=zero
          read(buffer,*,err=10) coor(1)
       else
          if(uniq == 0 .and. coor(1) >= zero) cycle
       end if
       if(coor(1) < zero) then
          iwork=coor(1)
          exit
       end if

       allocate(tmp_data,stat=status)
       if(status /= 0) call error_handler("QMMM1:read: failed TMP_DATA allocation")
       tmp_data%data%number=coor(1)
       tmp_data%data%x=coor(2)
       tmp_data%data%y=coor(3)
       tmp_data%data%z=coor(4)
       tmp_data%data%charge=zero
       tmp_data%data%i_uq=uniq
       nullify(tmp_data%next_data)
       current_data%next_data=>tmp_data
       current_data=>tmp_data

       n_at=n_at+1
    end do

    N_gx_atoms=n_at
    allocate(gx_qmmm(N_gx_atoms),stat=status)
    if(status /= 0) call error_handler("QMMM1: failed GX_QMMM allocation")

    current_data=>first_data
    do i=1,N_gx_atoms
       del=>current_data
       current_data=>current_data%next_data
       gx_qmmm(i)%number=current_data%data%number
       gx_qmmm(i)%x=current_data%data%x
       gx_qmmm(i)%y=current_data%data%y
       gx_qmmm(i)%z=current_data%data%z
       gx_qmmm(i)%charge=current_data%data%charge
       gx_qmmm(i)%i_uq=current_data%data%i_uq
       if(i > 1) deallocate(del)
    end do
    nullify(current_data)

    call returnclose_iounit(gxf)

    allocate(gx_qmmm_grad(N_gx_atoms),stat=status)
    if(status /= 0) call error_handler("QMMM1:read: failed GX_QMMM_GRAD allocation")
    do i=1,N_gx_atoms
       gx_qmmm_grad(i)%x=zero
       gx_qmmm_grad(i)%y=zero
       gx_qmmm_grad(i)%z=zero
    end do

    Energy_qmmm=zero; Energy_qm=zero; Energy_mm=zero

    !transfer of QM atoms from qmmm area to  unique_atoms area
    i=0
    do i_ua=1,N_unique_atoms
       counter_equal=1
       do
          if(counter_equal>unique_atoms(i_ua)%n_equal_atoms) exit
          i=i+1
          if(counter_equal==1) then
             unique_atoms(i_ua)%position_first_ea(1)=gx_qmmm(i)%x
             unique_atoms(i_ua)%position_first_ea(2)=gx_qmmm(i)%y
             unique_atoms(i_ua)%position_first_ea(3)=gx_qmmm(i)%z
          endif
          counter_equal=counter_equal+1
       enddo
    enddo! i_ua=1,N_unique_atoms
    N_gx_qm_atoms=i
    if (N_gx_qm_atoms >= N_gx_atoms) &
         call error_handler(" &
         & QMMM1: number of QM atoms >= total number of atoms")

    unique_atom_iwork=int(-iwork,i4_kind)

    return

10  call error_handler("QMMM1: failed GXFILE reading in")
  end subroutine read_gx_qmmm
  !*********************************************************

  !*********************************************************
  function read_qmmm1()
    !for QM+MM
    !------------ modules used -----------------------------------
    use species_module
    use slab_module
    !------------ Declaration of local variables -----------------
    logical :: read_qmmm1

    integer (i4_kind) :: i
    !------------ Executable code --------------------------------
    if(n_species /= N_gx_atoms) call error_handler( &
         "QMMM1: Number of atoms in molmech.inp does not coincide with "// &
         "those in gxfile")

    do i=1,n_species
       atoms_cart(i)%r(1)=gx_qmmm(i)%x*b2a
       atoms_cart(i)%r(2)=gx_qmmm(i)%y*b2a
       atoms_cart(i)%r(3)=gx_qmmm(i)%z*b2a
    end do

    if(lattice_calc) then
       call cart2frac(n_species,vect,atoms_cart,atoms_frac)
       call frac2cart(n_species,vect,atoms_frac,atoms_cart)
    else if(slab_calc) then
       do i=1,n_species
          call cart2frac_slab(atoms_cart(i)%r,atoms_frac(i)%r)
          call frac2cart_slab(atoms_frac(i)%r,atoms_cart(i)%r)
       end do
    end if

    read_qmmm1 = .true.
  end function read_qmmm1
  !******************************************************************

  !*********************************************************
  subroutine qmmm2pc()
    !------------ modules used -----------------------------------
    use datatype, only: pointcharge_type
    use pointcharge_module
    use operations_module, only : operations_solvation_effect
    use solv_cavity_module, only : with_pc,fixed_pc
    use species_module, only : n_species,atoms_cart,atoms
    use qmmm_interface_module, only : qm_mm_1
    !------------ Declaration of local variables -----------------
    type gx2pc
       type(pointcharge_type) :: data
       type(gx2pc), pointer :: next_data
    end type gx2pc
    type(gx2pc), target :: first_data
    type(gx2pc), pointer :: current_data, tmp_data, del

    integer :: i,j,n,i_uq_pc,i_uq_pc_old,n_eq
    integer :: status
    !------------ Executable code --------------------------------

    do i=1,n_species
       j=atoms_cart(i)%type
       gx_qmmm(i)%charge=atoms(j)%charge
       gx_qmmm(i)%name=atoms(j)%name
    end do

    current_data=>first_data
    nullify(current_data%next_data)

    n=0; i_uq_pc_old=0
    do i=N_gx_qm_atoms+1,N_gx_atoms
       i_uq_pc=gx_qmmm(i)%i_uq
       if(i_uq_pc == i_uq_pc_old) then
          n_eq=n_eq+1
          tmp_data%data%N_equal_charges=n_eq
          cycle
       else
          allocate(tmp_data,stat=status)
          if(status /= 0) call error_handler("QMMM1:pc: failed TMP_DATA allocation")
          n=n+1
          n_eq=1
          tmp_data%data%position_first_ec(1)=gx_qmmm(i)%x
          tmp_data%data%position_first_ec(2)=gx_qmmm(i)%y
          tmp_data%data%position_first_ec(3)=gx_qmmm(i)%z
          tmp_data%data%Z=gx_qmmm(i)%charge
          tmp_data%data%name=gx_qmmm(i)%name
          tmp_data%data%N_equal_charges=n_eq
          nullify(tmp_data%next_data)
          current_data%next_data=>tmp_data
          current_data=>tmp_data
          i_uq_pc_old=i_uq_pc
       end if
    end do

    if(pointcharge_N /= 0) call error_handler( &
         "QMMM1: Pointcharges were defined twice")
    pointcharge_N=n
    allocate( pointcharge_array(pointcharge_N), stat=status )
    if (status .ne. 0) call error_handler( &
         "QMMM1: pointcharge_array allocate failed.")
!!$    if((.not.options_spin_orbit).and.pointcharge_N>0) pseudopot_present=.true.

    current_data=>first_data
    do i=1,pointcharge_N
       del=>current_data
       current_data=>current_data%next_data
       pointcharge_array(i)%position_first_ec=current_data%data%position_first_ec
       pointcharge_array(i)%N_equal_charges=current_data%data%N_equal_charges
       pointcharge_array(i)%name=current_data%data%name
       pointcharge_array(i)%Z=current_data%data%Z
       pointcharge_array(i)%C=zero
       if(i > 1) deallocate(del)
    end do
    nullify(current_data)

    if(operations_solvation_effect) then
       with_pc=.true.
       fixed_pc=qm_mm_1
    end if

  end subroutine qmmm2pc
  !*********************************************************

  !*********************************************************
  subroutine QMfield_at_mm_points()
    !
    ! Executed by all workers. Called from main_master().
    !
    use pointcharge_module
    use elec_static_field_module
    use density_data_module, only: density_data_free
    use options_module, only: options_integrals_on_file
    use integralstore_module, only: integralstore_deallocate_pcm
    implicit none
    ! *** end of interface ***

    real(r8_kind), allocatable :: P(:,:)
    integer :: status,n,i,j,k


    n=N_gx_atoms-N_gx_qm_atoms
    allocate(P(3,n),stat=status)
    if(status /= 0) call error_handler("QMMM1:efield: failed P allocation")

    k=0
    do i=1,pointcharge_N
       do j=1,pointcharge_array(i)%N_equal_charges
          k=k+1
          P(:,k)=pointcharge_array(i)%position(:,j)
       end do
    end do

    calc_normal = .false.
    call fill_surf_points(P,N)

    deallocate(P,stat=status)
    if(status /= 0) call error_handler("QMMM1:efield: failed P deallocation")

    if(N_surface_points /= pointcharge_N) &
         call error_handler("QMMM1: wrong transfer from pointcharge array to surface_points")

    call surf_points_grad_information()
    call field_calculate ()
    call start_read_field_e()
    call get_field_nuc()

    ! This procedure runs in a parallel context. The following cleanup
    ! is executed by all workers:
    call density_data_free()
    if (.not. options_integrals_on_file()) then
       call integralstore_deallocate_pcm()
    end if
  end subroutine QMfield_at_mm_points
  !*********************************************************

  !*********************************************************
  subroutine qm_grads_to_qmmm1()
    !
    ! This sub  is called from a  parallel context but  on master only
    ! (see main_gradient()). It cannot use communication as the slaves
    ! are not doing anything "during  this time"! Check this next time
    ! it is running.
    !
    use unique_atom_module
    use energy_calc_module, only : get_energy
    use gradient_data_module, only : gradient_cartesian
    use elec_static_field_module, only : deallocate_field,E_ele,E_nuc,N_surface_points
    use elec_static_field_module, only : bounds_free_field,destroy_field_file
    use elec_static_field_module, only : surf_points_gradinfo_dealloc,dealloc_surf_points
    use solv_electrostat_module, only : F_solv_pc,dealloc_solv_pc
    use operations_module, only : operations_solvation_effect
    use qmmm_interface_module, only : qm_mm_1
    implicit none
    ! *** end of interface ***

    real(r8_kind) :: energy,pc_charge
    integer :: i,j,k,l
    character(len=2) :: qmmm_type

    ABORT("check!")
    call get_energy(tot=energy)
    Energy_qm=energy

    k=0
    do i=1,N_unique_atoms
       l=unique_atoms(i)%moving_atom
       do j=1,unique_atoms(i)%n_equal_atoms
          k=k+1
          if(l==0) then
             gx_qmmm_grad(k)%x = 0.0_r8_kind
             gx_qmmm_grad(k)%y = 0.0_r8_kind
             gx_qmmm_grad(i)%z = 0.0_r8_kind
          else
             gx_qmmm_grad(k)%x = gradient_cartesian(i)%m(1,j)
             gx_qmmm_grad(k)%y = gradient_cartesian(i)%m(2,j)
             gx_qmmm_grad(k)%z = gradient_cartesian(i)%m(3,j)
          endif
       end do
    end do

    if(k /= N_gx_qm_atoms) call error_handler("QMMM1: something wrong :-(")

    if( .not. qm_mm_1) then
       do i=1,N_surface_points
          l=size(E_ele(i)%m,2)
          do j=1,l
             k=k+1
             pc_charge=gx_qmmm(k)%charge
             gx_qmmm_grad(k)%x = pc_charge*(E_ele(i)%m(1,j)+E_nuc(i)%m(1,j))
             gx_qmmm_grad(k)%y = pc_charge*(E_ele(i)%m(2,j)+E_nuc(i)%m(2,j))
             gx_qmmm_grad(k)%z = pc_charge*(E_ele(i)%m(3,j)+E_nuc(i)%m(3,j))
             if(operations_solvation_effect) then
                gx_qmmm_grad(k)%x = gx_qmmm_grad(k)%x+F_solv_pc(i)%m(1,j)
                gx_qmmm_grad(k)%y = gx_qmmm_grad(k)%y+F_solv_pc(i)%m(2,j)
                gx_qmmm_grad(k)%z = gx_qmmm_grad(k)%z+F_solv_pc(i)%m(3,j)
             end if
          end do
       end do

       if(operations_solvation_effect) call dealloc_solv_pc
       call surf_points_gradinfo_dealloc()
       call bounds_free_field()
       call destroy_field_file()
       call dealloc_surf_points()
       call deallocate_field()
    end if

    !Output QM energy and gradients
    write(output_unit,*) ""
    write(output_unit,*) ""
    write(output_unit,*) "******************** QM+MM ***********************"
    write(output_unit,'(a21,f25.15)') "Energy after QM run: ",energy_qm
    write(output_unit,*) ""
    write(output_unit,'(a23)') "Gradients after QM run:"
    do i=1,N_gx_atoms
       if(i <= N_gx_qm_atoms) then
          qmmm_type="QM"
       else
          qmmm_type="MM"
       end if
       write(output_unit,'(a3,3f25.15,2x,a2)') &
            gx_qmmm(i)%name,gx_qmmm_grad(i)%x,gx_qmmm_grad(i)%y,gx_qmmm_grad(i)%z,qmmm_type
    end do
    write(output_unit,*) "**************************************************"
  end subroutine qm_grads_to_qmmm1
  !*********************************************************

  !*********************************************************
  subroutine mm_grads_to_qmmm1()
    !------------ modules used -----------------------------------
    use energy_and_forces_module, only: E_total,Grad
    use qmmm_interface_module, only : qm_mm_1
    !------------ Declaration of local variables -----------------
    integer :: i
    character(len=2) :: qmmm_type
    real(r8_kind) :: G(3)
    !------------ Executable code --------------------------------

    Energy_mm=E_total/h2kJm

    write(output_unit,*) "******************** QM+MM ***********************"
    write(output_unit,'(a21,f25.15)') "Energy after MM run: ",energy_mm
    write(output_unit,*) ""
    write(output_unit,'(a23)') "Gradients after MM run:"
    do i=1,N_gx_atoms
       G=zero
       if(i <= N_gx_qm_atoms) then
          qmmm_type="QM"
          G=Grad(:,i)
       else
          qmmm_type="MM"
          if(.not. qm_mm_1) G=Grad(:,i)
       end if
       write(output_unit,'(a3,3f25.15,2x,a2)') gx_qmmm(i)%name, &
            G(1)/(h2kJm*a2b),G(2)/(h2kJm*a2b),G(3)/(h2kJm*a2b),qmmm_type
    end do
    write(output_unit,*) "**************************************************"
    Energy_qmmm=Energy_qm+Energy_mm

    do i=1,N_gx_atoms
       G=zero
       if(i <= N_gx_qm_atoms) then
          G=Grad(:,i)
       else
          if(.not. qm_mm_1) G=Grad(:,i)
       end if
       gx_qmmm_grad(i)%x=gx_qmmm_grad(i)%x+G(1)/(h2kJm*a2b)
       gx_qmmm_grad(i)%y=gx_qmmm_grad(i)%y+G(2)/(h2kJm*a2b)
       gx_qmmm_grad(i)%z=gx_qmmm_grad(i)%z+G(3)/(h2kJm*a2b)
    end do

    write(output_unit,*) "******************** QM+MM ***********************"
    write(output_unit,'(a19,f25.15)') "Total QMMM Energy: ",energy_qmmm
    write(output_unit,*) ""
    write(output_unit,'(a21)') "Total QMMM Gradients:"
    do i=1,N_gx_atoms
       if(i <= N_gx_qm_atoms) then
          qmmm_type="QM"
       else
          qmmm_type="MM"
       end if
       write(output_unit,'(a3,3f25.15,2x,a2)') &
            gx_qmmm(i)%name,gx_qmmm_grad(i)%x,gx_qmmm_grad(i)%y,gx_qmmm_grad(i)%z,qmmm_type
    end do
    write(output_unit,*) "**************************************************"

  end subroutine mm_grads_to_qmmm1
  !*********************************************************

  !*********************************************************
  subroutine write_gx_qmmm()
    !------------ modules used -----------------------------------
    !------------ Declaration of local variables -----------------
    integer :: gxf,ios,status
    integer, parameter :: max_gx_atoms=1000
    character(len=256) :: buffer
    integer :: i,j,N_at
    real(r8_kind), allocatable :: z_num(:),xyz(:,:)
    integer, allocatable :: uq_num(:),at_num(:),zmat(:,:),numx(:,:)
    !------------ Executable code --------------------------------

    allocate(z_num(max_gx_atoms),xyz(3,max_gx_atoms),uq_num(max_gx_atoms), &
         at_num(max_gx_atoms),zmat(3,max_gx_atoms),numx(3,max_gx_atoms), &
         stat=status)
    if(status /= 0) call error_handler("QMMM1: failed allocation of GX dimentions")

    gxf = openget_iounit(file=trim(inpfile('gxfile')), status='unknown', &
         form='formatted')

    i=0
    do
       i=i+1
       read(gxf,'(a256)') buffer
!!$       read(buffer,'(f7.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4)',iostat=ios) &
       read(buffer,*,iostat=ios) &
            z_num(i),xyz(1:3,i),uq_num(i),at_num(i),zmat(1:3,i),numx(1:3,i)
!!$       if(ios /= 0) read(buffer,'(f7.2)') z_num(i)
       if(ios /= 0) read(buffer,*) z_num(i)
       if(z_num(i) < 0) exit
    end do
    N_at=i-1

    rewind gxf

    j=0
    do i=1,N_at
       if(uq_num(i) /= 0) then
          j=j+1
          xyz(1,i)=gx_qmmm(j)%x
          xyz(2,i)=gx_qmmm(j)%y
          xyz(3,i)=gx_qmmm(j)%z
       end if
       write(gxf,'(f5.2,3(2x,f21.12),2i4,2x,3I4,2X,3I4)') &
            z_num(i),xyz(1:3,i),uq_num(i),at_num(i),zmat(1:3,i),numx(1:3,i)
    end do
    write(gxf,'(f5.0,3(2x,f21.12),2i4,2x,3I4,2X,3I4)') z_num(N_at+1),0.0,0.0,0.0,0,0,0,0,0,0,0,0

    write(gxf,'(2F24.12,2x,3I5)') energy_qmmm,energy_qmmm,0,0,0
    do i=1,N_gx_atoms
       write(gxf,'(I5,5x,3F17.12)') i,gx_qmmm_grad(i)%x,gx_qmmm_grad(i)%y,gx_qmmm_grad(i)%z
    end do

    call returnclose_iounit(gxf)

    deallocate(z_num,xyz,uq_num,at_num,zmat,numx,stat=status)
    if(status /= 0) call error_handler("QMMM1: failed allocation of GX dimentions")

    deallocate(gx_qmmm,gx_qmmm_grad,stat=status)
    if(status /= 0) call error_handler("QMMM1: failed allocation of GX_QMMM and GX_QMMM_GRAD")

  end subroutine write_gx_qmmm
end module qmmm1_interface_module

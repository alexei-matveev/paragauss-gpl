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
module qmmm_interface_module
#include "def.h"
  !------------ Modules used --------------------------------------
  use type_module
  use iounitadmin_module
  use filename_module
  use string_qmmm_module

  implicit none
  private
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  !------------ General qmmm interface variables
  logical, public :: efp     !see efp_module
  logical, public :: imomm
  logical, public :: qm_mm   !see qmmm1_interface_module
  logical, public :: qm_mm_1 !see qmmm1_interface_module
  integer, public, parameter :: mm_run=0, qmmm_read_input=-1
  integer, public, parameter :: imomm_mm_small=2, imomm_mm_large=3
  integer, public, parameter :: qm_mm_run=4
  integer, public :: qm_mm_1_task
  !--------------------------------------------------------------

  integer(i4_kind), public, parameter :: max_atom = 1000
  integer(i4_kind), public :: n_total_at
  integer(i4_kind), public :: n_qm_at

  type,public :: gx_data
     real(r8_kind)          ::  charge
     real(r8_kind)          ::  x,y,z
     integer(i4_kind)       ::  i_unique,i_atom
     integer(i4_kind)       ::  zmat1,zmat2,zmat3
     integer(i4_kind)       ::  num1,num2,num3
  end type gx_data
  type(gx_data),public :: gx(max_atom)
  type(gx_data),public :: gx_qm(max_atom)

  type,public :: gradient
     real(r8_kind) :: x,y,z
  end type gradient
  type(gradient),public :: grad_qm(max_atom)
  type(gradient),public :: grad_mm2(max_atom)
  type(gradient),public :: grad_mm3(max_atom)

  real(r8_kind),public :: energy_qm, energy_mm2, energy_mm3
  !------------ public functions and subroutines ------------------
  public qmmm_read,qmmm_write,read_qmmm_input,sum_up_grads_and_write_gx
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  character(len=10) :: method
  character(len=10) :: df_method="imomm" ! "qm+mm", "qm+mm_pc"
  namelist /qmmm/ method

  real(r8_kind) :: energy_total

  type(gradient) :: grad_total(max_atom)
  type(gradient) :: grad_local(max_atom)
  type(gradient) :: grad_link(max_atom)

  integer(i4_kind), parameter :: string_length = 120
  integer(i4_kind), parameter :: max_link_atom = 20   ! Maximum no. of link atoms in the QM/MM system
  integer(i4_kind), parameter :: max_nml_opt = 100    ! number of input lines of namelist variable
                                                      ! that attached below gx file, i.e. when
                                                      ! using free valence coordinates
  type nml_link
     integer                ::  qm_atom
     integer                ::  mm_atom
     logical                ::  fixed_distance
     real(r8_kind)          ::  distance
     real(r8_kind)          ::  ratio
     integer                ::  atom_type       ! MM atom type for link atom.
  end type nml_link
  type (nml_link) :: linkinfo(max_link_atom)

  integer(i4_kind) :: n_atom
  integer(i4_kind) :: n_atom_gx
  integer(i4_kind) :: n_qm_atom
  integer(i4_kind) :: n_mm_atom
  integer(i4_kind) :: n_nml_input_opt

  character(len=string_length) :: nml_input_opt(max_nml_opt) ! Namelist input for optimizer
  real(r8_kind) :: energy

  integer(i4_kind) :: n_link_atoms
  integer(i4_kind) :: qm_atom
  integer(i4_kind) :: mm_atom
  logical :: fixed_distance
  real(r8_kind) :: distance
  real(r8_kind) ::  ratio
  integer(i4_kind) ::  atom_type

  integer(i4_kind) :: df_qm_atom=0
  integer(i4_kind) :: df_mm_atom=0
  logical :: df_fixed_distance=.true.
  real(r8_kind) :: df_distance=2.0786996_r8_kind  ! (au)= 1.10 Angstroem
  real(r8_kind) ::  df_ratio=0.708_r8_kind ! Morocuma famous values
  integer(i4_kind) ::  df_atom_type=5
  real(r8_kind) :: df_dummy_atom=99.00_r8_kind

  namelist /link_atom_number/  n_link_atoms
  namelist /link_atom/         qm_atom,mm_atom,fixed_distance,distance,ratio,atom_type
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine qmmm_read()
    !------------ Modules used ------------------------------------
    use input_module
    use inp_out_module, only : upcase, check_string
    use operations_module, only : operations_qm_mm_new
    !------------ Declaration of formal parameters ----------------
    !== End of interface ==========================================
    !------------ Declaration of local variables ------------------
    integer :: unit, status
    !------------ Executable code ---------------------------------

    method=df_method
    imomm=.false.
    qm_mm=.false.
    qm_mm_1=.false.
    efp=.false.

    if(input_line_is_namelist("qmmm")) then
       unit = input_intermediate_unit()
       call input_read_to_intermediate
       read(unit, nml=qmmm, iostat=status)
       if (status .gt. 0) call input_error("qmmm_read: namelist qmmm")
    endif

    if(.not. operations_qm_mm_new) then
       imomm=.false.
       qm_mm=.false.
       qm_mm_1=.false.
       efp=.false.
       return
    end if

    call upcase(method)

    efp=check_string(method,"EFP")
    if(efp) then
#ifndef WITH_EFP
       if(operations_qm_mm_new) call error_handler("EFP method was not compiled. Recompile -DWITH_EFP")
#endif
       return
    end if
    imomm=(index(method,"imomm") /= 0 .or. index(method,"IMOMM") /= 0)
    if(imomm) return
    qm_mm_1=(index(method,"qm+mm_pc") /= 0 .or. index(method,"QM+MM_PC") /= 0)
    if(qm_mm_1) return
    qm_mm=(index(method,"qm+mm") /= 0 .or. index(method,"QM+MM") /= 0)
    if(qm_mm) return

    call error_handler("qmmm_read: what type of QMMM calculations do you want to run?")

  end subroutine qmmm_read
  !****************************************************************

  !****************************************************************
  subroutine qmmm_write(iounit)
    !------------ Modules used ------------------------------------
    use echo_input_module
    use operations_module, only: operations_echo_input_level
    !------------ Declaration of formal parameters ----------------
    integer(kind=i4_kind), intent(in) :: iounit
    !== End of interface ==========================================
    !------------ Declaration of local variables ------------------
    !------------ Executable code ---------------------------------

    word_format = '("    ",a," = ",a15 :" # ",a)'

    call start("QMMM","QMMM_WRITE",iounit,operations_echo_input_level)
    call word("METHOD",method,df_method)
    call stop()

  end subroutine qmmm_write
  !****************************************************************

  !****************************************************************
  subroutine init_arrs()

    integer(i4_kind) :: counter

    do counter = 1,max_atom
       gx(counter)%x=0.0_r8_kind
       gx(counter)%y=0.0_r8_kind
       gx(counter)%z=0.0_r8_kind
       gx_qm(counter)%x=0.0_r8_kind
       gx_qm(counter)%y=0.0_r8_kind
       gx_qm(counter)%z=0.0_r8_kind
       grad_total(counter)%x=0.0_r8_kind
       grad_total(counter)%y=0.0_r8_kind
       grad_total(counter)%z=0.0_r8_kind
       grad_local(counter)%x=0.0_r8_kind
       grad_local(counter)%y=0.0_r8_kind
       grad_local(counter)%z=0.0_r8_kind
       grad_qm(counter)%x=0.0_r8_kind
       grad_qm(counter)%y=0.0_r8_kind
       grad_qm(counter)%z=0.0_r8_kind
       grad_mm2(counter)%x=0.0_r8_kind
       grad_mm2(counter)%y=0.0_r8_kind
       grad_mm2(counter)%z=0.0_r8_kind
       grad_mm3(counter)%x=0.0_r8_kind
       grad_mm3(counter)%y=0.0_r8_kind
       grad_mm3(counter)%z=0.0_r8_kind
    end do

  end subroutine init_arrs
  !****************************************************************

  !****************************************************************
  subroutine read_qmmm_input(saved_opt_step) !!PUBLIC

    real(r8_kind),intent(out) :: saved_opt_step

    character(len=*),parameter :: master_input   = 'master.input'
    character(len=*),parameter :: gxfile_all     = 'gxfile' !'gx'
!!$    character(len=*),parameter :: qmmm_all       = 'qmmm.all'
    logical :: exist_flag


    call init_arrs()

    call read_master_input(master_input)

    n_atom_gx = n_atom
    n_total_at=n_atom

    ! Check whether the gxfile.all exists or not.
    ! If this file exist, read it in.
    ! If this file does not exist, then append -1.0 at the end of gx information.
    inquire(exist=exist_flag, file=trim(inpfile(gxfile_all)))
    if (exist_flag) then
       call read_gx(gxfile_all)
       saved_opt_step      =  gx(n_atom_gx+1)%charge
    else
       gx(n_atom+1)%charge = -1.0
       saved_opt_step      = -1.0
    end if

    ! Write global intermediate file "qmmm.all".
!!$    call write_qmmm_all(qmmm_all)

    gx_qm(1:n_qm_atom)=gx(1:n_qm_atom)

    call data_preparing()

    n_qm_at=n_qm_atom+n_link_atoms

  end subroutine read_qmmm_input
  !****************************************************************

  !****************************************************************
  subroutine read_master_input(input_name)

    !This subroutine reads master input and return many global variables
    !that can determine further operations. Be careful !

    ! NOTE: This subroutine is independent and can be used any time.

    integer(i4_kind) :: i_unit
    integer(i4_kind) :: i,j,next
    character(len=*),intent(in)    :: input_name
    character(len=*),parameter     :: qm_pattern     = '%QM_REGION'
    character(len=*),parameter     :: int_pattern    = '%INTERFACE'
    character(len=*),parameter     :: enviro_pattern = '%ENVIRONMENT'
    character(len=*),parameter     :: link_pattern   = '%LINK_ATOMS'
    character(len=*),parameter     :: opt_pattern    = '%OPTIMIZER'
    character(len=string_length)   :: record,string
    integer(i4_kind) :: n_check

    i_unit=openget_iounit(trim(inpfile(input_name)), &
         form='formatted', status='unknown')

    ! Set total number of atoms, number of QM atoms and MM atoms to zero
    n_atom    = 0
    n_qm_atom = 0
    n_mm_atom = 0

    ! Check qm_pattern, if found exit from the loop, otherwise stop.

    check_qm_pattern: do while(.true.)
       read(unit=i_unit,fmt='(a)',end=10,err = 90) record

       if (record(1:len_trim(qm_pattern)) .eq. qm_pattern ) then
          exit check_qm_pattern
       else
          cycle check_qm_pattern
       end if
10     call error_handler("QMMM_INTERFACE: "//trim(qm_pattern)//" not found")
90     call error_handler("QMMM_INTERFACE: an error reading the master.input file")
    end do check_qm_pattern

    ! Read QM input and check for interface pattern (int_pattern). If
    ! found then proceeds, otherwise stops.

    check_int_pattern: do while(.true.)
       read(unit=i_unit,fmt='(a)',end=20,err=100) record

       if (record(1:len_trim(int_pattern)) .eq. int_pattern) then
          exit check_int_pattern
       else
          n_atom    = n_atom + 1
          n_qm_atom = n_qm_atom + 1
          read(record,*) gx(n_atom)%charge, &
               gx(n_atom)%x,gx(n_atom)%y,gx(n_atom)%z, &
               gx(n_atom)%i_unique,gx(n_atom)%i_atom, &
               gx(n_atom)%zmat1,gx(n_atom)%zmat2,gx(n_atom)%zmat3, &
               gx(n_atom)%num1,gx(n_atom)%num2,gx(n_atom)%num3

          cycle check_int_pattern
       end if
20     call error_handler("QMMM_INTERFACE: "//trim(int_pattern)//" not found")
100    call error_handler("QMMM_INTERFACE: an error reading qm information of the master.input file")
    end do check_int_pattern

    ! Check whether QM input is empty data or not.
    if (n_qm_atom .eq.0 ) call error_handler("QMMM_INTERFACE: There is no QM atoms")

    ! Now read in interface atoms and include them into QM atom lists.
    ! Check for enviro_pattern while reading.

    check_enviro_pattern: do while(.true.)
       read(unit=i_unit,fmt='(a)',end=30,err = 110) record

       if (record(1:len_trim(enviro_pattern)) .eq. enviro_pattern) then
          exit check_enviro_pattern
       else
          n_atom    = n_atom + 1
          n_mm_atom = n_mm_atom + 1
          read(record,*) gx(n_atom)%charge, &
               gx(n_atom)%x,gx(n_atom)%y,gx(n_atom)%z, &
               gx(n_atom)%i_unique,gx(n_atom)%i_atom, &
               gx(n_atom)%zmat1,gx(n_atom)%zmat2,gx(n_atom)%zmat3, &
               gx(n_atom)%num1,gx(n_atom)%num2,gx(n_atom)%num3
          cycle check_enviro_pattern
       end if
30     call error_handler("QMMM_INTERFACE: "//trim(enviro_pattern)//" not found")
110    call error_handler( &
            "QMMM_INTERFACE: an error reading interface atom information of the master.input file")
    end do check_enviro_pattern

    ! Check whether Interface input is empty data or not.
    if (n_mm_atom .eq. 0 ) call error_handler("QMMM_INTERFACE: There is no Interface atoms")

    ! Save number of interface atom here, so that it can be checked
    ! later whether MM input data is empty or not.
    n_check = n_mm_atom

    ! Now read in MM atoms and check for link_pattern.
    check_link_pattern: do while(.true.)
       read(unit=i_unit,fmt='(a)',end=40,err = 120) record

       if (record(1:len_trim(link_pattern)) .eq. link_pattern) then
          exit check_link_pattern
       else
          n_atom    = n_atom +1
          n_mm_atom = n_mm_atom + 1
          read(record,*) gx(n_atom)%charge, &
               gx(n_atom)%x,gx(n_atom)%y,gx(n_atom)%z, &
               gx(n_atom)%i_unique,gx(n_atom)%i_atom, &
               gx(n_atom)%zmat1,gx(n_atom)%zmat2,gx(n_atom)%zmat3, &
               gx(n_atom)%num1,gx(n_atom)%num2,gx(n_atom)%num3
          cycle check_link_pattern
       end if
40     call error_handler("QMMM_INTERFACE: "//trim(link_pattern)//" not found")
120    call error_handler(&
            "QMMM_INTERFACE: an error reading environment &atom information of the master.input file")
    end do check_link_pattern

    ! Check whether MM input is empty data or not.
    if (n_mm_atom .eq. n_check) call error_handler("QMMM_INTERFACE: There is no environment atoms")

    ! Now read the namelists of link atoms information
    n_link_atoms = 0                   ! Set zero the number of link atom.
    read(unit=i_unit,nml=link_atom_number)
    if (n_link_atoms .eq. 0 ) call error_handler("QMMM_INTERFACE: There is no link atoms")

    ! After knowing how many link atoms, now read in the link information.

    do i = 1,n_link_atoms
       ! Set default to the namelists variable first.
       qm_atom        = df_qm_atom
       mm_atom        = df_mm_atom
       fixed_distance = df_fixed_distance
       distance       = df_distance
       ratio          = df_ratio
       atom_type      = df_atom_type

       ! Read link information.
       read(unit=i_unit,nml=link_atom,err=130,end=130)

       ! Check some errors in link informations
       if (qm_atom .eq. 0 .or. mm_atom .eq. 0 ) call error_handler( &
            "QMMM_INTERFACE: error in the Namelist &LINK_ATOM. QM or MM number can not be 0")
       if (qm_atom .ge. mm_atom) call error_handler( &
            "QMMM_INTERFACE: error in the Namelist &LINK_ATOM. QM number can not greater than or equal to MM number")

       ! Save link information in saved array.
       linkinfo(i)%qm_atom        = qm_atom
       linkinfo(i)%mm_atom        = mm_atom
       linkinfo(i)%fixed_distance = fixed_distance
       linkinfo(i)%distance       = distance
       linkinfo(i)%ratio          = ratio
       linkinfo(i)%atom_type      = atom_type
    end do

    ! Now find opt_pattern.

    check_opt_pattern: do while(.true.)
       read(unit=i_unit,fmt='(a)',end=60, err = 140) record

       if (record(1:len_trim(opt_pattern)) .eq. opt_pattern) then
          exit check_opt_pattern
       else
          cycle check_opt_pattern
       end if
140    call error_handler("QMMM_INTERFACE: Cant read past link info")
    end do check_opt_pattern

    ! After opt_pattern was found, now store the input lines.
    j = 0
    read_opt_input: do i = 1,max_nml_opt
       read(unit=i_unit,fmt='(a)',end=60,err = 150) record

       ! Check whether it is empty line or not. We do not store empty lines.
       next = 1
       call gettext(record,string,next)
       if (next .eq. 1) then
          cycle read_opt_input
       else
          j = j + 1
          nml_input_opt(j) = record
       end if
    end do read_opt_input
    n_nml_input_opt = j      ! Number of Namelist lines for Optimizer

60  call returnclose_iounit(i_unit)
    return

130 call error_handler("QMMM_INTERFACE: Unable to read link atom information")
150 call error_handler("QMMM_INTERFACE: Cant read optimizer info")

  end subroutine read_master_input
  !****************************************************************

  !****************************************************************
  subroutine read_gx(gx_name)

    ! This subroutine returns values of array "gx" and also total number
    ! of atoms as "n_atom_gx"

    character(len=*), intent(in)  :: gx_name  !gx file name
    integer(i4_kind) :: gx_unit,i
    character(len=string_length) :: record

    gx_unit=openget_iounit(trim(inpfile(gx_name)), &
         form='formatted', status='unknown')

    i=0

    loop1: do
       i = i+1
       read(gx_unit,'(a)',err=30) record     ! Read a line of input

       read(record,*) gx(i)%charge    ! Read the first colum of input from the input line

       ! Determine if the end of gx_data has been reached or not
       if (gx(i)%charge .le. 0 ) then
          n_atom_gx = i - 1
          exit loop1
       else
          read(record,*) gx(i)%charge, &
               gx(i)%x,gx(i)%y,gx(i)%z, &
               gx(i)%i_unique,gx(i)%i_atom, &
               gx(i)%zmat1,gx(i)%zmat2,gx(i)%zmat3, &
               gx(i)%num1,gx(i)%num2,gx(i)%num3
       end if
    end do loop1

    call returnclose_iounit(gx_unit)

    return

30  call error_handler( &
         "QMMM_INTERFACE: an error during reading "//trim(gx_name)//" file")

  end subroutine read_gx
  !****************************************************************

  !****************************************************************
  subroutine data_preparing()

    integer(i4_kind) :: link_atom_type   ! MM atom type for link atom
    integer(i4_kind) :: i
    integer(i4_kind) :: i_qm_atom, i_mm_atom
    real(r8_kind) :: r_x,r_y,r_z
    real(r8_kind) :: r
    real(r8_kind) :: dist_ratio

    ! Determine number of QM atoms + Link atoms and calculate gx
    ! and connectivity information for link atoms.

    n_atom_gx = n_qm_atom

    do i = 1,n_link_atoms
       i_qm_atom = linkinfo(i)%qm_atom
       i_mm_atom = linkinfo(i)%mm_atom
       link_atom_type = linkinfo(i)%atom_type

       r_x = gx(i_mm_atom)%x - gx(i_qm_atom)%x
       r_y = gx(i_mm_atom)%y - gx(i_qm_atom)%y
       r_z = gx(i_mm_atom)%z - gx(i_qm_atom)%z

       r = sqrt( (r_x*r_x)+(r_y*r_y)+(r_z*r_z) )

       ! Determine distance ratio.
       if (linkinfo(i)%fixed_distance) then
          dist_ratio = linkinfo(i)%distance/r
       else
          dist_ratio = linkinfo(i)%ratio
       end if

       ! Add linkatom information into the gx information series.

       n_atom_gx = n_atom_gx + 1

       gx_qm(n_atom_gx)%x = gx(i_qm_atom)%x + (dist_ratio*r_x)
       gx_qm(n_atom_gx)%y = gx(i_qm_atom)%y + (dist_ratio*r_y)
       gx_qm(n_atom_gx)%z = gx(i_qm_atom)%z + (dist_ratio*r_z)

       gx_qm(n_atom_gx)%charge   = 1.00
       gx_qm(n_atom_gx)%i_atom   = n_atom_gx

    end do

  end subroutine data_preparing
  !****************************************************************

  !****************************************************************
  subroutine sum_up_grads_and_write_gx() !Public

    integer(i4_kind) :: i,j,i_start,i_end,k
    real(r8_kind) :: g     ! g = (R2-R1)/(R3-R1)
    real(r8_kind) :: r31   ! R3-R1
    real(r8_kind) :: r31x,r31y,r31z ! Components of R31 on X,Y,Z
    integer(i4_kind) :: i_qm_atom, i_mm_atom

    ! First read gradients for MM(R0,R1,R2,R4)
    k=0
    do i=1,n_total_at
       if(gx(i)%charge == df_dummy_atom) cycle
       k=k+1
       grad_local(i)%x = grad_mm3(k)%x
       grad_local(i)%y = grad_mm3(k)%y
       grad_local(i)%z = grad_mm3(k)%z
    end do
    energy_total=energy_mm3

    ! Second, read gradients for QM(R0,R1,R2) and
    ! Sum up the gradients for R0 and R1.
    k=0
    do i=1,n_qm_at-n_link_atoms !(or n_qm_atom)
       if(gx_qm(i)%charge == df_dummy_atom) cycle
       k=k+1
       grad_local(i)%x = grad_local(i)%x + grad_qm(k)%x
       grad_local(i)%y = grad_local(i)%y + grad_qm(k)%y
       grad_local(i)%z = grad_local(i)%z + grad_qm(k)%z
    end do
    energy_total=energy_total + energy_qm

    ! Store gradients for link atoms (R2)
    i_start=n_qm_atom+1
    i_end=n_qm_at
    j = 0
    do i=i_start,i_end
       j=j+1
       if(gx_qm(j+i_start)%charge == df_dummy_atom) cycle
       k=k+1
       grad_link(j)%x = grad_qm(k)%x
       grad_link(j)%y = grad_qm(k)%y
       grad_link(j)%z = grad_qm(k)%z
    end do

    ! Third, read gradients for MM(R0,R1,R2) and
    ! Sum up the gradients for R0 and R1.
    k=0
    do i=1,n_qm_at-n_link_atoms !(or n_qm_atom)
       if(gx_qm(i)%charge == df_dummy_atom) cycle
       k=k+1
       grad_local(i)%x = grad_local(i)%x - grad_mm2(k)%x
       grad_local(i)%y = grad_local(i)%y - grad_mm2(k)%y
       grad_local(i)%z = grad_local(i)%z - grad_mm2(k)%z
    end do
    energy_total=energy_total - energy_mm2

    ! Sum up gradients for link atoms
    j = 0
    do i=i_start,i_end
       j=j+1
       if(gx_qm(j+i_start)%charge == df_dummy_atom) cycle
       k=k+1
       grad_link(j)%x = grad_link(j)%x - grad_mm2(k)%x
       grad_link(j)%y = grad_link(j)%y - grad_mm2(k)%y
       grad_link(j)%z = grad_link(j)%z - grad_mm2(k)%z
    end do

    ! Transfer saved gradients to grad_total arrays.
    do i = 1,n_total_at
       grad_total(i)%x = grad_local(i)%x
       grad_total(i)%y = grad_local(i)%y
       grad_total(i)%z = grad_local(i)%z
    end do

    ! Fourth, transfer gradients from link atoms to R1 and R3.
    do i = 1,n_link_atoms
       i_qm_atom = linkinfo(i)%qm_atom
       i_mm_atom = linkinfo(i)%mm_atom

       r31x = gx(i_mm_atom)%x - gx(i_qm_atom)%x
       r31y = gx(i_mm_atom)%y - gx(i_qm_atom)%y
       r31z = gx(i_mm_atom)%z - gx(i_qm_atom)%z

       r31 = sqrt((r31x*r31x)+(r31y*r31y)+(r31z*r31z))

       ! Determine distance ratio.
       if (linkinfo(i)%fixed_distance) then
          ! ------ Option (b) ---------
          g = linkinfo(i)%distance / r31

          grad_total(i_qm_atom)%x = grad_total(i_qm_atom)%x + (grad_link(i)%x * deriv21_b(1,1))    &
               + (grad_link(i)%y * deriv21_b(2,1)) + (grad_link(i)%z * deriv21_b(3,1))

          grad_total(i_qm_atom)%y = grad_total(i_qm_atom)%y + (grad_link(i)%x * deriv21_b(1,2))    &
               + (grad_link(i)%y * deriv21_b(2,2)) + (grad_link(i)%z * deriv21_b(3,2))

          grad_total(i_qm_atom)%z = grad_total(i_qm_atom)%z + (grad_link(i)%x * deriv21_b(1,3))    &
               + (grad_link(i)%y * deriv21_b(2,3)) + (grad_link(i)%z * deriv21_b(3,3))

          grad_total(i_mm_atom)%x = grad_total(i_mm_atom)%x + (grad_link(i)%x * deriv23_b(1,1))    &
               + (grad_link(i)%y * deriv23_b(2,1)) + (grad_link(i)%z * deriv23_b(3,1))

          grad_total(i_mm_atom)%y = grad_total(i_mm_atom)%y + (grad_link(i)%x * deriv23_b(1,2))    &
               + (grad_link(i)%y * deriv23_b(2,2)) + (grad_link(i)%z * deriv23_b(3,2))

          grad_total(i_mm_atom)%z = grad_total(i_mm_atom)%z + (grad_link(i)%x * deriv23_b(1,3))    &
               + (grad_link(i)%y * deriv23_b(2,3)) + (grad_link(i)%z * deriv23_b(3,3))
       else
          ! ------ Option (a) ---------
          g = linkinfo(i)%ratio

          grad_total(i_qm_atom)%x = grad_total(i_qm_atom)%x + (grad_link(i)%x * deriv21_a(1,1))    &
               + (grad_link(i)%y * deriv21_a(2,1)) + (grad_link(i)%z * deriv21_a(3,1))

          grad_total(i_qm_atom)%y = grad_total(i_qm_atom)%y + (grad_link(i)%x * deriv21_a(1,2))    &
               + (grad_link(i)%y * deriv21_a(2,2)) + (grad_link(i)%z * deriv21_a(3,2))

          grad_total(i_qm_atom)%z = grad_total(i_qm_atom)%z + (grad_link(i)%x * deriv21_a(1,3))    &
               + (grad_link(i)%y * deriv21_a(2,3)) + (grad_link(i)%z * deriv21_a(3,3))

          grad_total(i_mm_atom)%x = grad_total(i_mm_atom)%x + (grad_link(i)%x * deriv23_a(1,1))    &
               + (grad_link(i)%y * deriv23_a(2,1)) + (grad_link(i)%z * deriv23_a(3,1))

          grad_total(i_mm_atom)%y = grad_total(i_mm_atom)%y + (grad_link(i)%x * deriv23_a(1,2))    &
               + (grad_link(i)%y * deriv23_a(2,2)) + (grad_link(i)%z * deriv23_a(3,2))

          grad_total(i_mm_atom)%z = grad_total(i_mm_atom)%z + (grad_link(i)%x * deriv23_a(1,3))    &
               + (grad_link(i)%y * deriv23_a(2,3)) + (grad_link(i)%z * deriv23_a(3,3))
       end if
    end do

    call write_gx()

  contains
    !============================================================
    function e_vec(i)
      ! This function calculate the scalar value of the unit vector
      ! R3-R1 along the axis 1,2,3 (X,Y,Z).
      real(r8_kind) :: e_vec
      integer(i4_kind), intent(in) :: i

      select case(i)
      case (1)
         e_vec = r31x / r31
      case (2)
         e_vec = r31y / r31
      case (3)
         e_vec = r31z / r31
      end select

    end function e_vec
    !============================================================
    function delta(i,j)
      ! This function return value of Kronecker delta
      real(r8_kind) :: delta
      integer(i4_kind), intent(in) :: i,j

      if (i.eq.j) then
         delta = 1.0
      else
         delta = 0.0
      end if

    end function delta
    !============================================================
    function deriv21_a(i,j)
      ! This function returns derivative of d(R2)/d(R1) for option (a)
      real(r8_kind) :: deriv21_a
      integer(i4_kind), intent(in) :: i,j

      deriv21_a = (1-g)*delta(i,j)

    end function deriv21_a
    !============================================================
    function deriv23_a(i,j)
      ! This function returns derivative of d(R2)/d(R3) for option (a)
      real(r8_kind) :: deriv23_a
      integer(i4_kind), intent(in) :: i,j

      deriv23_a = g*delta(i,j)

    end function deriv23_a
    !============================================================
    function deriv21_b(i,j)
      ! This function returns derivative of d(R2)/d(R1) for option (b)
      real(r8_kind) :: deriv21_b
      integer(i4_kind), intent(in) :: i,j

      deriv21_b = delta(i,j) + ( g*e_vec(i)*e_vec(j) ) - g*delta(i,j)

    end function deriv21_b
    !============================================================
    function deriv23_b(i,j)
      ! This function returns derivative of d(R2)/d(R3) for option (b)
      real(r8_kind) :: deriv23_b
      integer(i4_kind), intent(in) :: i,j

      deriv23_b = -( g*e_vec(i)*e_vec(j) ) + g*delta(i,j)

    end function deriv23_b

  end subroutine sum_up_grads_and_write_gx
  !****************************************************************

  !****************************************************************
  subroutine write_gx

    integer(i4_kind) :: out_unit,i
    character(len=*),parameter :: gx_format =   '(f5.2,3(2x,f21.12),2i4,2x,3i4,2x,3i4)'
    character(len=*),parameter :: gx_format_1 = '(f6.1,3(2X,f21.12),2i4,2x,3I4,2X,3I4)'
    character(len=*),parameter :: ener_format = '(2f21.12,2x,3i5)'
    character(len=*),parameter :: grad_format = '(i5,5x,3f17.12)'

!!$    out_unit=openget_iounit(trim(data_dir)//"/gxfile", &
!!$         form='formatted', status='unknown')
    out_unit=openget_iounit(trim(inpfile("gxfile")), &
         form='formatted', status='unknown')

    rewind(out_unit)

    ! Write gx informations
    do i = 1,n_atom
       write(unit=out_unit,fmt=gx_format,err = 20) gx(i)%charge, &
            gx(i)%x,gx(i)%y,gx(i)%z, &
            gx(i)%i_unique,gx(i)%i_atom, &
            gx(i)%zmat1,gx(i)%zmat2,gx(i)%zmat3, &
            gx(i)%num1,gx(i)%num2,gx(i)%num3

    end do

!!$    write(unit=out_unit,fmt='(f5.0)',err=20) gx(n_atom+1)%charge
    write(unit=out_unit,fmt=gx_format_1) gx(n_atom+1)%charge, &
         0.0_r8_kind,0.0_r8_kind,0.0_r8_kind,0,0,0,0,0,0,0,0

    ! Write Namelist information for Optimizer.
    do i = 1,n_nml_input_opt
       write(unit=out_unit,fmt='(a)',err=30) nml_input_opt(i)
    end do

    ! Write energy informations
    write(unit=out_unit,fmt=ener_format,err=40) energy_total,energy_total, &
         n_atom,n_atom,n_atom

    ! Write gradient informations.

    do i = 1,n_atom
       if (gx(i)%charge .ne. df_dummy_atom) then
          write(unit=out_unit,fmt=grad_format,err=50) i, grad_total(i)
       end if
    end do

    call returnclose_iounit(out_unit)

    return

20  call error_handler("QMMM_INTERFACE:write_gx_grad_opt:couldnt write gx information")
30  call error_handler( &
         "QMMM_INTERFACE:write_gx_grad_opt:couldnt write optimizer namelists information")
40  call error_handler("QMMM_INTERFACE:write_gx_grad_opt:couldnt write energy information")
50  call error_handler("QMMM_INTERFACE:write_gx_grad_opt:couldnt write gradient information")

  end subroutine write_gx
  !****************************************************************

  !****************************************************************
end module qmmm_interface_module

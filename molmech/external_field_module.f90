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
module external_field_module
  !------------ Modules used --------------------------------------
  use type_module
  use inp_out_module
  use common_data_module
  use species_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  logical, public :: ext_field
  logical, public :: ext_pc

  integer(i4_kind), public :: n_ext_pc

  type, public :: pc_coord
     integer(i4_kind) :: type
     real(r8_kind) :: r(3)
  end type pc_coord

  type, public :: pc_species
     character(len=len_name) :: name
     real(kind=r8_kind) :: charge
     real(kind=r8_kind) :: r_coval
  end type pc_species

  type(pc_coord), allocatable, public :: pc_cart(:)
  type(pc_species), allocatable, public :: p_charges(:)

  integer(i4_kind) :: n_ext_pc_types
  !------------ public functions and subroutines ------------------
  public read_ext_field_type, write_ext_field_type_to_output, read_ext_field, &
       pc_name2type, shutdown_ext_pc
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  character(len=6) :: type
  character(len=6), parameter :: df_type="NONE" !"PC"

  namelist /external_field_type/ type
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  function read_ext_field_type()
    logical :: read_ext_field_type

    integer(i4_kind) :: i
    logical :: read_tasks

    type=df_type

    read_ext_field_type = .false.

    call  go_to_first_input_line
    read_tasks=find_namelist("&EXTERNAL_FIELD_TYPE",i)
    if(.not.read_tasks) return
    read(input_device,nml=external_field_type, end=100, err=200)

    call upcase(type)

    select case (trim(type))
    case("NONE")
       ext_field=.false.
       ext_pc=.false.
       n_ext_pc=0
    case("PC")
       ext_field=.true.
       ext_pc=.true.
    case default
       call error_handler("MolMech: Wrong external force keyword")
    end select

    read_ext_field_type = .true.
    return

100 read_ext_field_type = .false.
    return
200 call input_nm_error(0,"EXTERNAL_FIELD_TYPE") 

  end function  read_ext_field_type
  !****************************************************************

  !****************************************************************
  subroutine write_ext_field_type_to_output()

    character(len=76) :: message

    write(output_device,'(80("*"))')

    if(ext_pc) then
       message='External field: point charge array'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if

    write(output_device,'(80("*"),/)')

  end subroutine write_ext_field_type_to_output
  !****************************************************************

  !****************************************************************
  subroutine read_ext_field()

    integer(i4_kind) :: i,j,status
    logical :: read_pc_coor
    character(len=80) :: string
    character(len=6) :: number
    character(len=len_name) :: pc_nm
    real(kind=r8_kind) :: coor(3),Q_pc
    integer(i4_kind) :: pc_type
    integer(kind=i4_kind), allocatable :: pc_type_tmp(:)
    character(len=len_name), allocatable :: pc_nm_tmp(:)
    real(kind=r8_kind), allocatable :: q_tmp(:)

    type inp_data
       type (pc_coord) :: data
       type (inp_data), pointer :: next_data
    end type inp_data
    type (inp_data), target :: first_data
    type (inp_data), pointer :: current_data, tmp_data, del  

    if(ext_pc) then
       read_pc_coor=find_word("&PC_COOR",i)
       if(.not.read_pc_coor) call error_handler &
            ("MolMech: No external PCs")

       allocate(pc_type_tmp(200),pc_nm_tmp(200),Q_tmp(200),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed PC_TYPE_TMP allocation")
       n_ext_pc_types=0
       
       current_data=>first_data
       nullify(current_data%next_data)

       i=0
       l1: do
          i=i+1
          read(input_device,'(a80)', err=100, end=200) string
          if(check_string(string,"&")) goto 200
          if(check_string(string,"/")) then
             call upcase(string)
             if(check_string(string,"/PC_COOR")) exit l1
             if(.not.check_string(string,"/PC_COOR")) goto 300
          end if
          read (string,*,err=400) pc_nm, coor(1:3), Q_pc

          pc_type=get_pc_type()
          allocate(tmp_data)
          tmp_data%data=pc_coord(pc_type,coor)
          nullify(tmp_data%next_data)
          current_data%next_data =>tmp_data
          current_data => tmp_data
       end do l1

       n_ext_pc=i-1
       allocate(pc_cart(n_ext_pc), stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed PC_CART allocation")

       current_data => first_data
       do j=1,n_ext_pc
          del=>current_data
          current_data=>current_data%next_data
          pc_cart(j)=current_data%data
          if(j > 1) deallocate(del)
       end do

       allocate(p_charges(n_ext_pc_types),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed P_CHARGES allocation")
       do j=1,n_ext_pc_types
          p_charges(j)=pc_species(pc_nm_tmp(j),q_tmp(j),r_ext_pc)
       end do

       deallocate(pc_type_tmp,pc_nm_tmp,q_tmp, stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed PC_TYPE_TMP allocation")

       return
       
100    write(number,'(i4)') i
       call error_handler("MolMech: Wrong "//trim(number)// &
            " line into PC_COOR block:"//achar(10)//trim(string))
200    call error_handler("MolMech: Not terminated PC_COOR block")
300    call error_handler("MolMech: Wrong definition of PC_COOR block")
400    write(number,'(i4)') i 
       call error_handler("MolMech: Coordinates of "//trim(number)// &
            " PC has been read in with error:"//achar(10)//trim(string))
    end if
    
  contains
    
    function get_pc_type()

      integer(kind=i4_kind) :: get_pc_type
      
      integer(kind=i4_kind) :: ig

      if(n_ext_pc_types==0) then
         n_ext_pc_types=1
         pc_nm_tmp(1)=pc_nm
         q_tmp(1)=Q_pc
         pc_type_tmp(1)=n_species_types+n_ext_pc_types
         get_pc_type=pc_type_tmp(1)
      else
         do ig=1,n_ext_pc_types
            if(pc_nm == pc_nm_tmp(ig)) then  
               if(Q_pc == Q_tmp(ig)) then
                  get_pc_type=pc_type_tmp(ig)
                  return
               else
                  write(number,'(i4)') i
                  call error_handler("MolMech: PC_COOR block: Check line "//trim(number)//"." &
                       //achar(10)//"Wrong correspondence between NAME and CHARGE")
               end if
            end if
         end do
         n_ext_pc_types=n_ext_pc_types+1
         pc_nm_tmp(n_ext_pc_types)=pc_nm
         q_tmp(n_ext_pc_types)=Q_pc
         pc_type_tmp(n_ext_pc_types)=n_species_types+n_ext_pc_types
         get_pc_type=pc_type_tmp(n_ext_pc_types)
      end if

    end function get_pc_type

  end subroutine read_ext_field
  !******************************************************************

  !******************************************************************
  function pc_name2type(pc_nm)

    integer(kind=i4_kind) :: pc_name2type
    character(len=len_name) :: pc_nm

    integer(kind=i4_kind) :: i

    pc_name2type=0_i4_kind
    do i=1,n_ext_pc_types
       if(trim(pc_nm) == trim(p_charges(i)%name)) then
          pc_name2type=i+n_species_types
          exit
       end if
    end do

  end function pc_name2type
  !******************************************************************

  !******************************************************************
  subroutine shutdown_ext_pc()

    integer(i4_kind) :: status
    
    deallocate(pc_cart,p_charges, stat=status)
    if(status /= 0) call error_handler( &
         "MolMech: failed PC_CART or. P_CHARGES deallocation")

  end subroutine shutdown_ext_pc

end module external_field_module

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
module potentials_module
  !------------ Modules used --------------------------------------
  use type_module
  use iounitadmin_module
  use filename_module
  use common_data_module
  use inp_out_module
  use species_module
  use tasks_main_options_module
  use external_field_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  type, public :: poten_data
     integer(kind=i4_kind) :: id        ! global identificator of the potential
     integer(kind=i4_kind) :: bond_type ! bonded or nonbonded
     integer(kind=i4_kind) :: Nb_type   ! N-body type
     integer(kind=i4_kind) :: n_atom    ! number atoms participated
     integer(kind=i4_kind) :: sp_type(10)  ! the list of species participating in the interaction
     integer(kind=i4_kind) :: n_param   ! number of parameters
     real(kind=r8_kind) :: param(n_parameter)    ! the list of potential parameters
  end type poten_data
  type(poten_data), allocatable, public :: poten(:)
  integer(kind=i4_kind), public :: n_in_pot

  type, public :: table
     integer(kind=i4_kind) :: n_points
     real(kind=r8_kind) :: step
     real(kind=r8_kind), pointer :: r(:), f(:), d2f(:), f1(:), d2f1(:)
  end type table
  type(table), public, allocatable :: user_poten(:),user_poten_tmp(:)
  integer(kind=i4_kind), public :: n_user_pot

  integer(kind=i4_kind), public :: n_2b,n_3b,n_4b,n_ss,n_bb,n_sb,n_st,n_2b_n,n_cs

  type, public :: list_of_potentials
     character(len=10) :: name(n_poten)   !name of potential cannot exceed 10 characters        
     integer(kind=i4_kind) :: type(n_poten)  ! 2 (2-body), 3 (3-body), 4 (4-body), 
     ! 5 (cross terms), 6 (core-shell)
     integer(kind=i4_kind) :: n_atoms(n_poten) ! number of atoms participated
     integer(kind=i4_kind) :: num_par(n_poten) ! number of parameters
     integer(kind=i4_kind) :: bonded_type(n_poten) ! 0-bonded, 1-nonbonded (for crossed terms
     !this parameter serves only to distinguish potentials each other)
  end type list_of_potentials
  type(list_of_potentials), public :: poten_list

  !------------ public functions and subroutines ------------------
  public potential_list_init, read_potential, write_potential_to_output, &  
       convert_ff_parameters, shutdown_poten_data
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  character(len=10) :: pot_name 
  character(len=len_name) :: atom_name(10) 
  real(kind=r8_kind) :: ff_parameter(n_parameter)
  character(len=200) :: read_data

  character(len=20) :: df_pot_name=' '
  character(len=len_name) :: df_atom_name(10)=(/'      ','      ','      ', &
       '      ','      ','      ','      ','      ','      ','      '/)
  real(kind=r8_kind) :: df_ff_parameter(n_parameter)=(/zero,zero,zero,zero,zero,zero,zero,zero, &
       zero,zero,zero,zero,zero,zero,zero/)
  character(len=200) :: df_read_data='NO'

  namelist /potential/ pot_name,atom_name,ff_parameter,read_data
  !------------ Subroutines ---------------------------------------
contains
  !****************************************************************
  subroutine potential_list_init()

    poten_list%name=(/'HARM_STR  ','MORZE     ','QUART_STR ','BUCK      ','L_J       ', &
         'HARM_BND  ','QUART_BND ','SIX_BND   ','HARM_TRS  ','TRIPL_COS ','STR_STR   ', &
         'BND_BND   ','STR_BND   ','STR_TRS   ','CORE_SHELL','BCK_SHORT '/)
    poten_list%type=       (/2,2,2,2,2,3,3,3,4,4,5,5,5,5,6,2/)
    poten_list%n_atoms=    (/2,2,2,2,2,3,3,3,4,4,3,5,3,4,2,3/)
    poten_list%num_par=    (/2,3,4,3,2,2,4,6,2,3,3,3,4,2,1,3/)
    poten_list%bonded_type=(/0,0,0,1,1,0,0,0,0,0,1,2,3,4,0,0/)

  end subroutine potential_list_init
  !****************************************************************

  !****************************************************************
  function read_potential()

    logical :: read_potential

    integer(kind=i4_kind) :: id,id1,type_a(10),nb,it,ib
    logical :: yes
    character(len=len_name) :: name(10)
    character(len=1) :: c_s(10)
    integer(kind=i4_kind) :: count,status,i,j,k,n_at,n_bt,n_t,n_np,ii
    integer(kind=i4_kind) :: k2,k3,k4,kss,kbb,ksb,kst,k2n,kcs
    real(kind=r8_kind) :: rconv

    type inp_data
       type (poten_data) :: data
       type (inp_data), pointer :: next_data
    end type inp_data
    type (inp_data), target :: first_data
    type (inp_data), pointer :: current_data, tmp_data, del

    if(trim(length_unit) == "ANGSTROM") rconv=one
    if(trim(length_unit) == "BOHR") rconv=b2a

    current_data=>first_data
    nullify(current_data%next_data)

    count=0_i4_kind
    n_2b=0; n_3b=0; n_4b=0; n_ss=0; n_bb=0; n_sb=0; n_st=0; n_2b_n=0; n_cs=0
    id=n_poten
    call go_to_first_input_line()
    l1:do
       pot_name=df_pot_name
       atom_name=df_atom_name
       ff_parameter=df_ff_parameter
       read_data=df_read_data

       read_potential=find_namelist("&POTENTIAL",ii)
       if(.not.read_potential) goto 100
       read(input_device,nml=potential, end=100, err=200)
       call upcase(pot_name)
       !definition global potential id
       if (check_string(pot_name,"USER")) then
          if(.not.allocated(user_poten_tmp)) then
             allocate(user_poten_tmp(100),stat=status)
             if(status /= 0) call error_handler( &
                  "MolMech: failed USER_POTEN_TMP allocation")
          end if
          id1=id+1
          id=id1
          if(trim(read_data)=="NO") goto 200
          call read_table(id1-n_poten,n_at,n_bt,n_t,n_np)
       else
          yes=.false.
          do i=1,n_poten
             if (check_string(pot_name,trim(poten_list%name(i)))) then
                id1=i
                n_at=poten_list%n_atoms(id1)
                n_bt=poten_list%bonded_type(id1)
                n_t=poten_list%type(id1)
                n_np=poten_list%num_par(id1)
                yes=.true.
                exit
             end if
          end do
          if (.not.yes) goto 200
       end if

       ! converting atom names into atom types
       i=0
       do 
          if(check_string(atom_name(i+1),"      ").or.i+1 > n_at) exit
          i=i+1
          if(check_string(atom_name(i)," S").or.check_string(atom_name(i)," s")) then
             c_s(i)="S"
          else
             c_s(i)="C"
          end if
          call name_without_cs(atom_name(i),name(i))
          type_a(i)=name2type(name(i),c_s(i))
          if(type_a(i)==0 .and. ext_field .and. ext_pc) type_a(i)=pc_name2type(name(i))
!!$          if(type_a(i)==0 .and. trim(poten_list%name(id1)) /= "BCK_SHORT") goto 200
          if(type_a(i) > n_species_types .and.  trim(poten_list%name(id1)) == "BCK_SHORT") &
               call error_handler("MolMech: BCK_SHORT potential still cannot be applied to"// &
               achar(10)//"regular species - external pc interaction")
          if(type_a(i) > n_species_types .and.  (id1 >= 9 .and. id1 <=15)) &
               call error_handler("MolMech: 4- and many-body potentials still cannot be applied to"// &
               achar(10)//"regular species - external pc interactions")
       end do
       nb=i
       if(nb < n_at) goto 200

       !ignore dummy potentials
       j=n_at
       do k=1,n_at 
          if(type_a(k) == 0) then
             if(k==3 .and. trim(poten_list%name(id1)) == "BCK_SHORT") exit
             j=0
             exit
          end if
       end do
       if(j == 0) cycle l1

       if(n_t==6) then
          i=type_a(1); j=type_a(2)
          if(i == j) goto 200
          if(get_c_s(atoms(i)%c_s) == "C" .and. get_c_s(atoms(j)%c_s) == "C") goto 200 
          if(get_c_s(atoms(i)%c_s) == "S" .and. get_c_s(atoms(j)%c_s) == "S") goto 200 
       end if

       !ignore potencials acting between external PC
       j=n_at
       do k=1,n_at 
          if(type_a(k) <= n_species_types) then
             j=0
             exit
          end if
       end do
       if(j == n_at) cycle l1

       
       if(n_t==2 .and. n_bt==0) n_2b=n_2b+1
       if(n_t==2 .and. n_bt==1) n_2b_n=n_2b_n+1
       if(n_t==3 .and. n_bt==0) n_3b=n_3b+1
       if(n_t==4 .and. n_bt==0) n_4b=n_4b+1
       if(n_t==5) then
          if(n_bt==1) n_ss=n_ss+1
          if(n_bt==2) n_bb=n_bb+1
          if(n_bt==3) n_sb=n_sb+1
          if(n_bt==4) n_st=n_st+1
       end if
       if(n_t==6) n_cs=n_cs+1

       allocate(tmp_data,stat=status)
       if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(poten)")
       tmp_data%data=poten_data(id1,n_bt,n_t,n_at,type_a,n_np,ff_parameter)
       nullify(tmp_data%next_data)
       current_data%next_data =>tmp_data
       current_data => tmp_data

       count=count+1
    end do l1
100 n_in_pot=count
    n_user_pot=id-n_poten


    if(n_in_pot == 0_i4_kind) then
       read_potential=.false.
    else
       read_potential=.true.

       allocate(poten(n_in_pot), stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed POTEN allocation")

       current_data => first_data

       k2=0
       k3=k2+n_2b
       k4=k3+n_3b
       kss=k4+n_4b; kbb=kss+n_ss; ksb=kbb+n_bb; kst=ksb+n_sb;
       k2n=kst+n_st
       kcs=k2n+n_2b_n
       do i=1,n_in_pot
          del=>current_data
          current_data => current_data%next_data
          it=current_data%data%Nb_type
          ib=current_data%data%bond_type
          if(it==2 .and. ib==0) then
             k2=k2+1
             poten(k2)=current_data%data
          else if(it==2 .and. ib==1) then
             k2n=k2n+1
             poten(k2n)=current_data%data
          else if(it==3) then
             k3=k3+1
             poten(k3)=current_data%data
          else if(it==4) then
             k4=k4+1
             poten(k4)=current_data%data
          else if(it==5) then
             if(ib==1) then
                kss=kss+1
                poten(kss)=current_data%data
             end if
             if(ib==2) then
                kbb=kbb+1
                poten(kbb)=current_data%data
             end if
             if(ib==3) then
                ksb=ksb+1
                poten(ksb)=current_data%data
             end if
             if(ib==4) then
                kst=kst+1
                poten(kst)=current_data%data
             end if
          else if(it==6) then
             kcs=kcs+1
             poten(kcs)=current_data%data
          end if
          if(i > 1) deallocate(del)
       end do
       nullify(current_data)

       if(n_user_pot > 0) then
          allocate(user_poten(n_user_pot),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed USER_POTEN allocation")
          do i=1,n_user_pot
             do k=1,n_in_pot
                if(n_poten+i == poten(k)%id) exit
             end do
             user_poten(i)%n_points=user_poten_tmp(i)%n_points
             user_poten(i)%step=user_poten_tmp(i)%step
             j=user_poten(i)%n_points
             allocate(user_poten(i)%r(j),user_poten(i)%f(j), &
                  user_poten(i)%d2f(j),stat=status)
             if(status /= 0) call error_handler("MolMech: Failed USER_POTEN%F allocation")
             if(poten(k)%nb_type == 2 .and. poten(k)%bond_type == 1) then
                allocate(user_poten(i)%f1(j),user_poten(i)%d2f1(j), stat=status)
                if(status /= 0) call error_handler("MolMech: Failed USER_POTEN%F1 allocation")
             end if
             user_poten(i)%r=user_poten_tmp(i)%r
             user_poten(i)%f=user_poten_tmp(i)%f
             user_poten(i)%d2f=user_poten_tmp(i)%d2f
             if(poten(k)%nb_type == 2 .and. poten(k)%bond_type == 1) then
                user_poten(i)%f1=user_poten_tmp(i)%f1
                user_poten(i)%d2f1=user_poten_tmp(i)%d2f1
             end if

             user_poten(i)%step=rconv*user_poten(i)%step
             user_poten(i)%r=rconv*user_poten(i)%r

             deallocate(user_poten_tmp(i)%r,user_poten_tmp(i)%f, &
                  user_poten_tmp(i)%d2f,stat=status)
             if(status /= 0) call error_handler("MolMech: Failed USER_POTEN_TMP%F deallocation")             
             if(poten(k)%nb_type == 2 .and. poten(k)%bond_type == 1) then
                deallocate(user_poten_tmp(i)%f1,user_poten_tmp(i)%d2f1, stat=status)
                if(status /= 0) call error_handler("MolMech: Failed USER_POTEN_TMP%F1 deallocation")
             end if
          end do
          deallocate(user_poten_tmp,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed USER_POTEN_TMP allocation")

       endif
    end if
    return

200 call input_nm_error(count+1,"POTENTIAL")

  contains

    subroutine read_table(icount,n_atoms,bonded_type,type,num_par)
      !read in tabulated potential function (kJ/mol and angstrom - only, while !!!!!!!!!)

      integer(kind=i4_kind), intent(in) :: icount
      integer(kind=i4_kind) :: n_atoms,bonded_type,type,num_par
      integer(kind=i4_kind) :: n_points
      real(kind=r8_kind) :: step
      integer(kind=i4_kind) :: i,file_id,ii
      character(len=200) :: line_buf,file_name,message
      logical :: yes,read_nml

      namelist /ident_data/ n_atoms,bonded_type,type,num_par,n_points,step

      !default parameters for vdW interactions - 
      !only interaction what user can define in current version   
      type=2; bonded_type=1; n_atoms=2; num_par=2   

      i=1
      do
         if (read_data(i:i) == " ") then
            i=i+1
         else
            exit
         end if
      end do
      line_buf=read_data(i:200)
      if (check_string(line_buf,"/")) then
         file_name=trim(line_buf)
      else
         file_name=trim(inpfile(line_buf))
      end if

      inquire(file=file_name,exist=yes)
      if (yes) then
         file_id=openget_iounit(file_name, form='formatted', status='old')
      else
         message="MolMech: File "//trim(file_name)//" is absent"
         call error_handler(message)
      end if

      read_nml=find_nml(file_id,"&IDENT_DATA",ii)
      if(.not.read_nml) goto 10
      read(file_id, nml=ident_data, end=10, err=20)
      user_poten_tmp(icount)%n_points=n_points
      user_poten_tmp(icount)%step=step

      allocate(user_poten_tmp(icount)%r(n_points),user_poten_tmp(icount)%f(n_points), &
           user_poten_tmp(icount)%d2f(n_points),stat=status)
      if(status /= 0) call error_handler("MolMech: Failed USER_POTEN_TMP%F allocation")
      if(type==2 .and. bonded_type==1) then
         allocate(user_poten_tmp(icount)%f1(n_points),user_poten_tmp(icount)%d2f1(n_points),stat=status)
         if(status /= 0) call error_handler("MolMech: Failed USER_POTEN_TMP%F1 allocation")
      end if

      do i=1,n_points
         if(type==2 .and. bonded_type==1) then
            read(file_id,*, end=20, err=20) &
                 user_poten_tmp(icount)%r(i),user_poten_tmp(icount)%f(i),user_poten_tmp(icount)%d2f(i), &
                 user_poten_tmp(icount)%f1(i),user_poten_tmp(icount)%d2f1(i)
         else
            read(file_id,*, end=20, err=20) &
                 user_poten_tmp(icount)%r(i),user_poten_tmp(icount)%f(i),user_poten_tmp(icount)%d2f(i)
         end if
         if(i > 1) then
            if(user_poten_tmp(icount)%r(i)-user_poten_tmp(icount)%r(i-1) > step+small .or. &
                 user_poten_tmp(icount)%r(i)-user_poten_tmp(icount)%r(i-1) < step-small ) goto 30
         end if
      end do
      call returnclose_iounit(file_id)
      return

10    message="MolMech: there is no namelist IDENT_DATA within file "//trim(file_name)
      call error_handler(message)
20    message="MolMech: error occered during reading in "//trim(file_name)
      call error_handler(message)
30    write(message,'("MolMech: Interval between ",i5," and ",i5," points is not equal to STEP")') i,i-1
      call error_handler(message)

    end subroutine read_table

  end function read_potential
  !----------------------------------------------------------------
  function find_nml(f_device,word,i)

    integer(kind=i4_kind), intent(in) :: f_device
    character(len=*), intent(in) :: word    
    integer(kind=i4_kind), intent(out) :: i

    logical :: find_nml

    character(len=80) :: string, message
    character(len=4) :: number

    i=0
    do 
       i=i+1
       read(f_device,'(a80)', end=100, err=200) string
       call upcase(string)
       if (check_string(string,word)) then
          find_nml= .true.
          exit
       end if
    end do
    backspace f_device
    return

100 find_nml= .false.
    return
200 write(number,'(i4)') i
    message = trim(" MolMech: The error reading in line "//trim(number)// &
         " of table file")
    call error_handler(message)

  end function find_nml
  !****************************************************************

  !****************************************************************
  subroutine convert_ff_parameters()

    integer(kind=i4_kind) :: i
    real(kind=r8_kind) :: econv,rconv,aconv

    econv=one; rconv=one; aconv=deg2rad

    if(trim(energy_unit) == "KCAL/MOL") then 
       econv=c2J
    elseif(trim(energy_unit) == "EV") then
       econv=eV2kJm
    end if
    if(trim(length_unit) == "BOHR") rconv=b2a
    if(trim(angle_unit) == "RADIAN") aconv=one
    do i=1,n_in_pot
       select case(poten(i)%id)
       case(1)  !harm_str
          poten(i)%param(1)=econv*poten(i)%param(1)/rconv**2
          poten(i)%param(2)=rconv*poten(i)%param(2)
       case(2) !morze
          poten(i)%param(1)=econv*poten(i)%param(1)
          poten(i)%param(2)=poten(i)%param(2)/rconv
          poten(i)%param(3)=rconv*poten(i)%param(3)
       case(3) !quart_str
          poten(i)%param(1)=econv*poten(i)%param(1)/rconv**2
          poten(i)%param(2)=econv*poten(i)%param(2)/rconv**3
          poten(i)%param(3)=econv*poten(i)%param(3)/rconv**4
          poten(i)%param(4)=rconv*poten(i)%param(4)
       case(4) !buck
          poten(i)%param(1)=econv*poten(i)%param(1)
          poten(i)%param(2)=rconv*poten(i)%param(2)
          poten(i)%param(3)=econv*poten(i)%param(3)*rconv**6
       case(5) !l_j
          poten(i)%param(1)=econv*poten(i)%param(1)*rconv**12
          poten(i)%param(2)=econv*poten(i)%param(2)*rconv**6
       case(6) ! harm_bnd
          poten(i)%param(1)=econv*poten(i)%param(1)/aconv**2
       case(7) ! quart_bnd
          poten(i)%param(1)=econv*poten(i)%param(1)/aconv**2
          poten(i)%param(2)=econv*poten(i)%param(2)/aconv**3
          poten(i)%param(3)=econv*poten(i)%param(3)/aconv**4
       case(8) ! six_bnd
          poten(i)%param(1)=econv*poten(i)%param(1)/aconv**2
          poten(i)%param(2)=econv*poten(i)%param(2)/aconv**3
          poten(i)%param(3)=econv*poten(i)%param(3)/aconv**4
          poten(i)%param(4)=econv*poten(i)%param(4)/aconv**5
          poten(i)%param(5)=econv*poten(i)%param(5)/aconv**6
       case(9) ! harm_trs
          poten(i)%param(1)=econv*poten(i)%param(1)/aconv**2
       case(10) ! tripl_cos
          poten(i)%param(1:3)=econv*poten(i)%param(1:3)
       case(11) ! str_str
          poten(i)%param(1)=econv*poten(i)%param(1)/rconv**2
          poten(i)%param(2)=rconv*poten(i)%param(3)
          poten(i)%param(3)=rconv*poten(i)%param(2)
       case(12) ! bnd_bnd
          poten(i)%param(1)=econv*poten(i)%param(1)
       case(13) !str_bnd
          poten(i)%param(1)=econv*poten(i)%param(1)/rconv
          poten(i)%param(2:3)=rconv*poten(i)%param(2:3)
       case(14) ! str_trs
          poten(i)%param(1)=econv*poten(i)%param(1)/rconv
          poten(i)%param(2)=rconv*poten(i)%param(2)
       case(15) !core_shell
          poten(i)%param(1)=econv*poten(i)%param(1)/rconv**2
       case(16) !bck_short
          poten(i)%param(1)=econv*poten(i)%param(1)
          poten(i)%param(2)=rconv*poten(i)%param(2)
          poten(i)%param(3)=econv*poten(i)%param(3)*rconv**6
       end select
    end do

  end subroutine convert_ff_parameters
  !****************************************************************

  !****************************************************************
  subroutine write_potential_to_output()

    integer(kind=i4_kind) :: i,j,k,l,na,np
    character(len=len_name) :: at_nm(10)
    character(len=1) :: c_s(10)
    character(len=12) ::  pot_nm, bnd_ty

    write(output_device,'(80("-"))')

    do i=1,n_in_pot
       if(i==1.and.n_2b/=0) write(output_device,'("Stratching",/,"==========")')
       if(i==n_2b+1.and.n_3b/=0) write(output_device,'("Bending",/,"=======")')
       if(i==n_2b+n_3b+1.and.n_4b/=0) write(output_device,'("Torsion",/,"=======")')
       if(i==n_2b+n_3b+n_4b+1.and.n_ss/=0) &
            write(output_device,'("Stratching-Stratching",/,"=====================")')
       if(i==n_2b+n_3b+n_4b+n_ss+1.and.n_bb/=0) &
            write(output_device,'("Bending_Bending",/,"===============")')
       if(i==n_2b+n_3b+n_4b+n_ss+n_bb+1.and.n_sb/=0) &
            write(output_device,'("Stratching-Bending",/,"==================")')
       if(i==n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+1.and.n_st/=0) &
            write(output_device,'("Stratching-Torsion",/,"==================")')
       if(i==n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+n_st+1) &
            write(output_device,'("Two-body nonbonding",/,"====================")')
       j=poten(i)%id
       if(j <= n_poten) then
          pot_nm=poten_list%name(j)
       else
          pot_nm="USER_DEF"
       end if
       if(poten(i)%bond_type == 1 .and. poten(i)%Nb_type == 2) then
          bnd_ty="non bonded"
       else
          bnd_ty="bonded"
       end if
       na=poten(i)%n_atom
       do k=1,na
          l=poten(i)%sp_type(k)
          if (l /= 0) then
             if (l <= n_species_types) then
                at_nm(k)=atoms(l)%name
                c_s(k)=get_c_s(atoms(l)%c_s)
                if(c_s(k) == "C") c_s(k)=" "
                if(c_s(k) == "S") c_s(k)="s"
             else
                at_nm(k)=p_charges(l-n_species_types)%name
                c_s(k)=" "
             end if
          else
             at_nm(k)=" "
             c_s(k)=" "
          end if
          at_nm(k)=trim(at_nm(k))//' '//trim(c_s(k))
       end do
       np=poten(i)%n_param
       write (output_device,'(1x,a14,1x,a14)') pot_nm, bnd_ty
!!$       write (output_device,'(10(1x,a6))') at_nm(1:na)
       write (output_device,'(10(1x,a6))') at_nm(1:na)
       write (output_device,'(15(1x,f15.6))') poten(i)%param(1:np)
       if(i/=n_in_pot) write(output_device,'(80("."),/)')
    end do
    write(output_device,'(80("-"),/)')

  end subroutine write_potential_to_output
  !****************************************************************

  !****************************************************************
  subroutine shutdown_poten_data()

    integer(i4_kind) :: status,i

    if(allocated(poten)) then
       deallocate(poten,stat=status)
       if(status /= 0) call error_handler("MolMech: failed POTEN deallocation")
    end if

    if(allocated(user_poten)) then
       do i=1,n_user_pot
          deallocate(user_poten(i)%r,user_poten(i)%f,user_poten(i)%d2f, &
               user_poten(i)%f1,user_poten(i)%d2f1, stat=status)
          if(status /= 0) call error_handler("MolMech: failed USER_POTEN data deallocation")
       end do
       deallocate(user_poten, stat=status)
       if(status /= 0) call error_handler("MolMech: failed USER_POTEN deallocation")
    end if

  end subroutine shutdown_poten_data
  !****************************************************************

  !****************************************************************
end module potentials_module






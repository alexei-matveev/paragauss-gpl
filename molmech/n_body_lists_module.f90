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
module n_body_lists_module
  !------------ Modules used -----------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use tasks_main_options_module
  use slab_module
  use species_module
  use potentials_module
  use external_field_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of public constants and variables -----
  logical, public :: nb_lists_from_input
  integer(kind=i4_kind), public, allocatable :: N_total(:) !define the total possible 
  !number of partners for Coulomb and van der Waals interactions 

  type, public :: lists
     integer(kind=i4_kind), pointer :: list(:,:)
  end type lists
  type(lists), public, allocatable :: two_body(:)
  type(lists), public, allocatable :: three_body(:)
  type(lists), public, allocatable :: four_body(:)
  type(lists), public, allocatable :: str_str(:)
  type(lists), public, allocatable :: bnd_bnd(:)
  type(lists), public, allocatable :: str_bnd(:)
  type(lists), public, allocatable :: str_trs(:)
  type(lists), public, allocatable :: core_shell(:)
  type(lists), public, allocatable :: bonded_spcs(:) !for vdw
  type(lists), public, allocatable :: bonded_spcs_c(:) !for ewald

  type, public :: list_nb
     integer(kind=i4_kind), pointer :: list(:,:)
     logical, pointer :: first_image(:)
  end type list_nb
  type(list_nb), public, allocatable :: non_bonded(:)
  type(list_nb), public, allocatable :: coulombic(:)

  type, public :: image_list
     integer(i4_kind), pointer :: list(:,:)
  end type image_list
  type(image_list), public, allocatable :: add_list(:) !for vdw
  type(image_list), public, allocatable :: add_list_c(:) !for ewald

  logical, public :: automatic_nb_lists
  integer(i4_kind), public  :: update_nb_list
  logical, public :: exclude_1_2
  logical, public :: ff3_style
  logical, public :: exclude_1_3
  logical, public :: exclude_1_4

  !------------ public functions and subroutines ---------------------
  public read_nb_options, read_Nb_lists, resorting_Nb_lists, autobuilding_nb_lists, &
       build_2b_nonbonded_lists, shutdown_interacting_lists, write_nb_options, &
       send_receive_n_total
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of private constants and variables ----
  logical :: df_automatic_nb_lists=.false.
  integer(i4_kind) :: df_update_nb_list=1
  logical :: df_exclude_1_2=.true.
  logical :: df_ff3_style=.false.
  logical :: df_exclude_1_3=.true.
  logical :: df_exclude_1_4=.false.

  namelist /nb_lists_options/ automatic_nb_lists,update_nb_list,exclude_1_2, ff3_style, &
       exclude_1_3, exclude_1_4
  !------------ Subroutines -----------------------------------------
contains
  !******************************************************************
  function read_nb_options()

    logical :: read_nb_options
    integer(i4_kind) :: i

    automatic_nb_lists=df_automatic_nb_lists
    update_nb_list=df_update_nb_list
    exclude_1_2=df_exclude_1_2
    ff3_style=df_ff3_style
    exclude_1_3=df_exclude_1_3
    exclude_1_4=df_exclude_1_4
    
    call  go_to_first_input_line
    read_nb_options=find_namelist("&NB_LISTS_OPTIONS",i)
    if(.not.read_nb_options) return
    read(input_device,nml=nb_lists_options, end=100, err=200)

    if(update_nb_list < 1 ) goto 200

    read_nb_options=.true.
    return

100 read_nb_options=.false.
    return
200 call input_nm_error(0,"NB_LISTS_OPTIONS")

  end function read_nb_options
  !******************************************************************

  !******************************************************************
  subroutine write_nb_options()

    logical :: yes
    character(len=76) :: message

    yes=.false.

    if(automatic_nb_lists) then
       yes=.true.
       message='Automatic generation of N-body lists. Only for 2-,3-,4-bodies and core-shell'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if(.not.exclude_1_2) then
       yes=.true.
       message='1-2 exclusion - false'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if(.not.exclude_1_3) then
       yes=.true.
       message='1-3 exclusion - false'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if(exclude_1_4) then
       yes=.true.
       message='1-4 exclusion - true'
       write(output_device,'(a2,a76,a2)') '* ',message,' *'
    end if
    if(yes) then
       write(output_device,'(80("*"),/)')
    else
       write(output_device,'(/)')
    end if

  end subroutine write_nb_options
  !******************************************************************

  !******************************************************************
  function read_Nb_lists()
    !read in N-body lists of bonded interactions from input

    logical :: read_Nb_lists

    integer(kind=i4_kind) :: i0,i1,n_2b_l
    integer(kind=i4_kind) :: n_3b_l,n_4b_l,n_ss_l,n_bb_l,n_sb_l,n_st_l,n_cs_l
    integer(kind=i4_kind) :: ty(10),num(10),num1(10)
    integer(kind=i4_kind), allocatable :: n_pot(:)
    logical :: yes(2)
    integer(kind=i4_kind) :: i,j,k,status,nl
    character(len=100) :: message

    type buf_data
       integer(kind=i4_kind) :: ind(10)
       integer(kind=i4_kind) :: n_pt
       type (buf_data), pointer :: next_buf
    end type buf_data
    type (buf_data), target :: first_buf
    type (buf_data), pointer :: current_buf, tmp_buf

    read_Nb_lists =.false.

    if(n_cs /= 0) then
!!$print*,"CORE_SHELL_LIST"
       yes(1)=find_word("&CORE_SHELL_LIST",i0)
       yes(2)=find_word("/CORE_SHELL_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 50
       if(i0 >= i1) goto 50
       n_cs_l=i1-i0-1
       yes(1)=find_word("&CORE_SHELL_LIST",i0)
       
       allocate(core_shell(n_cs),n_pot(n_cs),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed CORE_SHELL allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_cs_l: do i=1,n_cs_l
          read(input_device, *,err=100) num(1:2)

          if(num(1) == num(2)) then
             write(message,'(a42,2i5)') "MolMech: CORE_SHELL_LIST: wrong definition: ",num(1:2)
             call error_handler(trim(message))
          end if

          ty(1)=atoms_cart(num(1))%type
          ty(2)=atoms_cart(num(2))%type

          j_n_cs: do j=1,n_cs
             k=n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+n_st+n_2b_n+j

             yes(1)= ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2)
             yes(2)= ty(1)==poten(k)%sp_type(2) .and. ty(2)==poten(k)%sp_type(1)
             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_cs

             n_pot(j)=n_pot(j)+1
             allocate(tmp_buf)
             tmp_buf%ind(1)=num(1); tmp_buf%ind(2)=num(2)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf

         end do j_n_cs
       end do i_n_cs_l

       nl=0
       do i=1,n_cs
          nl=nl+n_pot(i)
          allocate(core_shell(i)%list(2,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                core_shell(i)%list(1,j)=current_buf%ind(1)
                core_shell(i)%list(2,j)=current_buf%ind(2)
!!$print*,'!',i,core_shell(i)%list(1:2,j)
             end if
          end do
       end do

       if(to_epe) then
          allocate(cs_pair(2,nl),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed CS_PAIR allocation")
          k=0
          do i=1,n_cs
             do j=1,n_pot(i)
                k=k+1
                cs_pair(1,k)=core_shell(i)%list(1,j)
                cs_pair(2,k)=core_shell(i)%list(2,j)
             end do
          end do
       end if

       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation -1")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    if(n_2b /= 0) then
!!$print*,"TWO_BODY_LIST"
       yes(1)=find_word("&TWO_BODY_LIST",i0)
       yes(2)=find_word("/TWO_BODY_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 100
       if(i0 >= i1) goto 100
       n_2b_l=i1-i0-1
       yes(1)=find_word("&TWO_BODY_LIST",i0)
       
       allocate(two_body(n_2b),n_pot(n_2b),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed TWO_BODY allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_2b_l: do i=1,n_2b_l
          read(input_device, *,err=100) num(1:2)

          if(num(1) == num(2)) then
             write(message,'(a42,2i5)') "MolMech: TWO_BODY_LIST: wrong definition: ",num(1:2)
             call error_handler(trim(message))
          end if

          if(num(1) <= n_species) ty(1)=atoms_cart(num(1))%type
          if(num(2) <= n_species) ty(2)=atoms_cart(num(2))%type
          if(num(1) > n_species) ty(1)=pc_cart(num(1)-n_species)%type
          if(num(2) > n_species) ty(2)=pc_cart(num(2)-n_species)%type

          j_n_2b: do j=1,n_2b

             yes(1)= ty(1)==poten(j)%sp_type(1) .and. ty(2)==poten(j)%sp_type(2)
             yes(2)= ty(1)==poten(j)%sp_type(2) .and. ty(2)==poten(j)%sp_type(1)
             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_2b

             if(ty(1) <= n_species_types) then 
                if(atoms(ty(1))%redfac /= one) atoms_cart(num(1))%neigh_sp=num(2)
             end if
             if(ty(2) <= n_species_types) then 
                if(atoms(ty(2))%redfac /= one) atoms_cart(num(2))%neigh_sp=num(1)
             end if

             n_pot(j)=n_pot(j)+1
             allocate(tmp_buf)
             tmp_buf%ind(1)=num(1); tmp_buf%ind(2)=num(2)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf
             
          end do j_n_2b
       end do i_n_2b_l

       do i=1,n_2b
          allocate(two_body(i)%list(2,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                two_body(i)%list(1,j)=current_buf%ind(1)
                two_body(i)%list(2,j)=current_buf%ind(2)
!!$print*,'!',i,two_body(i)%list(1:2,j)
             end if
          end do
       end do

       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation 0")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    if(n_3b /= 0) then
!!$print*,"THREE_BODY_LIST"
       yes(1)=find_word("&THREE_BODY_LIST",i0)
       yes(2)=find_word("/THREE_BODY_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 200
       if(i0 >= i1) goto 200
       n_3b_l=i1-i0-1
       yes(1)=find_word("&THREE_BODY_LIST",i0)
       
       allocate(three_body(n_3b),n_pot(n_3b),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed THREE_BODY allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_3b_l: do i=1,n_3b_l
          read(input_device,*,err=200) num(1:3)

          if(num(1) == num(2) .or. &
             num(1) == num(3) .or. &
             num(2) == num(3)) then
             write(message,'(a44,3i5)') "MolMech: THREE_BODY_LIST: wrong definition: ",num(1:3)
             call error_handler(trim(message))
          end if

          if(num(1) <= n_species) ty(1)=atoms_cart(num(1))%type
          if(num(2) <= n_species) ty(2)=atoms_cart(num(2))%type
          if(num(3) <= n_species) ty(3)=atoms_cart(num(3))%type
          if(num(1) > n_species) ty(1)=atoms_cart(num(1)-n_species)%type
          if(num(2) > n_species) ty(2)=atoms_cart(num(2)-n_species)%type
          if(num(3) > n_species) ty(3)=atoms_cart(num(3)-n_species)%type

          j_n_3b: do j=1,n_3b
             k=n_2b+j

             yes(1)= ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(3)
             yes(2)= ty(1)==poten(k)%sp_type(3) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(1)
             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_3b

             n_pot(j)=n_pot(j)+1
             allocate(tmp_buf)
             tmp_buf%ind(1)=num(1); tmp_buf%ind(2)=num(2); tmp_buf%ind(3)=num(3)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf
             
          end do j_n_3b
       end do i_n_3b_l

       do i=1,n_3b
          allocate(three_body(i)%list(3,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                three_body(i)%list(1,j)=current_buf%ind(1)
                three_body(i)%list(2,j)=current_buf%ind(2)
                three_body(i)%list(3,j)=current_buf%ind(3)
!!$print*,'!',i,three_body(i)%list(1:3,j)
             end if
          end do
       end do
       
       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation 1")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    if(n_4b /= 0) then
!!$print*,"FOUR_BODY_LIST"
       yes(1)=find_word("&FOUR_BODY_LIST",i0)
       yes(2)=find_word("/FOUR_BODY_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 300
       if(i0 >= i1) goto 300
       n_4b_l=i1-i0-1
       yes(1)=find_word("&FOUR_BODY_LIST",i0)

       allocate(four_body(n_4b),n_pot(n_4b),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed FOUR_BODY allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_4b_l: do i=1,n_4b_l
          read(input_device,*,err=300) num(1:4)

          if(num(1) == num(2) .or. &
             num(1) == num(3) .or. &
             num(1) == num(4) .or. &
             num(2) == num(3) .or. &
             num(2) == num(4) .or. &
             num(3) == num(4)) then
             write(message,'(a43,4i5)') "MolMech: FOUR_BODY_LIST: wrong definition: ",num(1:4)
             call error_handler(trim(message))
          end if

          ty(1:4)=atoms_cart(num(1:4))%type

          j_n_4b: do j=1,n_4b
             k=n_2b+n_3b+j

             yes(1)= ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(4)
             yes(2)= ty(1)==poten(k)%sp_type(4) .and. ty(2)==poten(k)%sp_type(3) .and. &
                  ty(3)==poten(k)%sp_type(2) .and. ty(4)==poten(k)%sp_type(1)
             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_4b

             if(yes(1)) then
                num1(1:4)=num(1:4)
             else if(yes(2)) then
                num1(1)=num(4)
                num1(2)=num(3)
                num1(3)=num(2)
                num1(4)=num(1)
             end if

             n_pot(j)=n_pot(j)+1

             allocate(tmp_buf)
             tmp_buf%ind(1:4)=num1(1:4)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf

          end do j_n_4b
       end do i_n_4b_l

       do i=1,n_4b
          allocate(four_body(i)%list(4,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                four_body(i)%list(1:4,j)=current_buf%ind(1:4)
!!$print*,'!',i,four_body(i)%list(1:4,j)
             end if
          end do
       end do

       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation 0")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    if(n_ss /= 0) then
!!$print*,"STR_STR_LIST"
       yes(1)=find_word("&STR_STR_LIST",i0)
       yes(2)=find_word("/STR_STR_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 400
       if(i0 >= i1) goto 400
       n_ss_l=i1-i0-1
       yes(1)=find_word("&STR_STR_LIST",i0)

       allocate(str_str(n_ss),n_pot(n_ss),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed STR_STR allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_ss_l: do i=1,n_ss_l
          read(input_device,*,err=400) num(1:3)

          if(num(1) == num(2) .or. &
             num(1) == num(3) .or. &
             num(2) == num(3)) then
             write(message,'(a41,4i5)') "MolMech: STR_STR_LIST: wrong definition: ",num(1:3)
             call error_handler(trim(message))
          end if

          ty(1:3)=atoms_cart(num(1:3))%type

          j_n_ss: do j=1,n_ss
             k=n_2b+n_3b+n_4b+j

             yes(1)= ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(3)
             yes(2)= ty(1)==poten(k)%sp_type(3) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(1)
             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_ss

             if(yes(1)) then
                num1(1:3)=num(1:3)
             else if(yes(2)) then
                num1(1)=num(3)
                num1(2)=num(2)
                num1(3)=num(1)
             end if

             n_pot(j)=n_pot(j)+1
             allocate(tmp_buf)
             tmp_buf%ind(1:3)=num1(1:3)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf
             
          end do j_n_ss
       end do i_n_ss_l

       do i=1,n_ss
          allocate(str_str(i)%list(3,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                str_str(i)%list(1:3,j)=current_buf%ind(1:3)
!!$print*,'!',i,str_str(i)%list(1:3,j)
             end if
          end do
       end do

       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation 0")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    if(n_bb /= 0) then
!!$print*,"BND_BND_LIST"
       yes(1)=find_word("&BND_BND_LIST",i0)
       yes(2)=find_word("/BND_BND_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 500
       if(i0 >= i1) goto 500
       n_bb_l=i1-i0-1
       yes(1)=find_word("&BND_BND_LIST",i0)

       allocate(bnd_bnd(n_bb),n_pot(n_bb),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed BND_BND allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_bb_l: do i=1,n_bb_l
          read(input_device,*,err=500) num(1:5)

          if(num(1) == num(2) .or. &
             num(1) == num(3) .or. &
             num(2) == num(3) .or. &
             num(4) == num(5) .or. &
             num(4) == num(3) .or. &
             num(5) == num(3) .or. &
            (num(1) == num(4) .and. num(2) == num(5)) .or. &
            (num(1) == num(5) .and. num(2) == num(4))) then
             write(message,'(a41,5i5)') "MolMech: BND_BND_LIST: wrong definition: ",num(1:5)
             call error_handler(trim(message))
          end if

          ty(1:5)=atoms_cart(num(1:5))%type
          j_n_bb: do j=1,n_bb
             k=n_2b+n_3b+n_4b+n_ss+j

             yes(1)= (ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(4) .and. &
                  ty(5)==poten(k)%sp_type(5)) .or. &
                  (ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(5) .and. &
                  ty(5)==poten(k)%sp_type(4)) .or. &
                  (ty(1)==poten(k)%sp_type(2) .and. ty(2)==poten(k)%sp_type(1) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(4) .and. &
                  ty(5)==poten(k)%sp_type(5)) .or. &
                  (ty(1)==poten(k)%sp_type(2) .and. ty(2)==poten(k)%sp_type(1) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(5) .and. &
                  ty(5)==poten(k)%sp_type(4))

             yes(2)= (ty(1)==poten(k)%sp_type(4) .and. ty(2)==poten(k)%sp_type(5) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(1) .and. &
                  ty(5)==poten(k)%sp_type(2)) .or. &
                  (ty(1)==poten(k)%sp_type(4) .and. ty(2)==poten(k)%sp_type(5) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(2) .and. &
                  ty(5)==poten(k)%sp_type(1)) .or. &
                  (ty(1)==poten(k)%sp_type(5) .and. ty(2)==poten(k)%sp_type(4) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(1) .and. &
                  ty(5)==poten(k)%sp_type(2)) .or. &
                  (ty(1)==poten(k)%sp_type(5) .and. ty(2)==poten(k)%sp_type(4) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(2) .and. &
                  ty(5)==poten(k)%sp_type(1)) 

             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_bb

             if(yes(1)) then
                num1(1:5)=num(1:5)
             else if(yes(2)) then
                num1(1)=num(5)
                num1(2)=num(4)
                num1(3)=num(3)
                num1(4)=num(2)
                num1(5)=num(1)
             end if

             n_pot(j)=n_pot(j)+1
             allocate(tmp_buf)
             tmp_buf%ind(1:5)=num1(1:5)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf
             
          end do j_n_bb
       end do i_n_bb_l

       do i=1,n_bb
          allocate(bnd_bnd(i)%list(5,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                bnd_bnd(i)%list(1:5,j)=current_buf%ind(1:5)
!!$print*,'!',i,bnd_bnd(i)%list(1:5,j)
             end if
          end do
       end do

       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation 0")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    if(n_sb /= 0) then
!!$print*,"STR_BND_LIST"
       yes(1)=find_word("&STR_BND_LIST",i0)
       yes(2)=find_word("/STR_BND_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 600
       if(i0 >= i1) goto 600
       n_sb_l=i1-i0-1
       yes(1)=find_word("&STR_BND_LIST",i0)

       allocate(str_bnd(n_sb),n_pot(n_sb),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed STR_BND allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_sb_l: do i=1,n_sb_l
          read(input_device,*,err=600) num(1:3)

          if(num(1) == num(2) .or. &
             num(1) == num(3) .or. &
             num(2) == num(3)) then
             write(message,'(a41,4i5)') "MolMech: STR_BND_LIST: wrong definition: ",num(1:3)
             call error_handler(trim(message))
          end if

          ty(1:3)=atoms_cart(num(1:3))%type

          j_n_sb: do j=1,n_sb
             k=n_2b+n_3b+n_4b+n_ss+n_bb+j

             yes(1)= ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(3)
             yes(2)= ty(1)==poten(k)%sp_type(3) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(1)
             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_sb

             if(yes(1)) then
                num1(1:3)=num(1:3)
             else if(yes(2)) then
                num1(1)=num(3)
                num1(2)=num(2)
                num1(3)=num(1)
             end if

             n_pot(j)=n_pot(j)+1
             allocate(tmp_buf)
             tmp_buf%ind(1:3)=num1(1:3)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf
             
          end do j_n_sb
       end do i_n_sb_l

       do i=1,n_sb
          allocate(str_bnd(i)%list(3,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                str_bnd(i)%list(1:3,j)=current_buf%ind(1:3)
!!$print*,'!',i,str_bnd(i)%list(1:3,j)
             end if
          end do
       end do

       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation 0")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    if(n_st /= 0) then
!!$print*,"STR_TRS_LIST"
       yes(1)=find_word("&STR_TRS_LIST",i0)
       yes(2)=find_word("/STR_TRS_LIST",i1)
       if(.not.yes(1) .or. .not.yes(2)) goto 700
       if(i0 >= i1) goto 700
       n_st_l=i1-i0-1
       yes(1)=find_word("&STR_TRS_LIST",i0)

       allocate(str_trs(n_st),n_pot(n_st),stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed STR_TRS allocation 1")

       current_buf=>first_buf
       nullify(current_buf%next_buf)

       n_pot=0
       i_n_st_l: do i=1,n_st_l
          read(input_device,*,err=700) num(1:4)

          if(num(1) == num(2) .or. &
             num(1) == num(3) .or. &
             num(1) == num(4) .or. &
             num(2) == num(3) .or. &
             num(2) == num(4) .or. &
             num(3) == num(4)) then
             write(message,'(a41,4i5)') "MolMech: STR_TRS_LIST: wrong definition: ",num(1:4)
             call error_handler(trim(message))
          end if

          ty(1:4)=atoms_cart(num(1:4))%type

          j_n_st: do j=1,n_st
             k=n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+j

             yes(1)= ty(1)==poten(k)%sp_type(1) .and. ty(2)==poten(k)%sp_type(2) .and. &
                  ty(3)==poten(k)%sp_type(3) .and. ty(4)==poten(k)%sp_type(4)
             yes(2)= ty(1)==poten(k)%sp_type(4) .and. ty(2)==poten(k)%sp_type(3) .and. &
                  ty(3)==poten(k)%sp_type(2) .and. ty(4)==poten(k)%sp_type(1)
             if(.not.yes(1) .and. .not.yes(2)) cycle j_n_st

             n_pot(j)=n_pot(j)+1

             allocate(tmp_buf)
             tmp_buf%ind(1:4)=num(1:4)
             tmp_buf%n_pt=j
             nullify(tmp_buf%next_buf)
             current_buf%next_buf=>tmp_buf
             current_buf=>tmp_buf

          end do j_n_st
       end do i_n_st_l

       do i=1,n_st
          allocate(str_trs(i)%list(4,n_pot(i)))
          current_buf=>first_buf
          j=0
          do 
             if (.not. associated(current_buf%next_buf)) exit
             current_buf=>current_buf%next_buf
             if(current_buf%n_pt == i) then
                j=j+1
                str_trs(i)%list(1:4,j)=current_buf%ind(1:4)
!!$print*,'!',i,str_trs(i)%list(1:4,j)
             end if
          end do
       end do

       deallocate(n_pot,stat=status)
       if(status /= 0) call error_handler( &
         "MolMech: failed N_POT allocation 0")

       read_Nb_lists=.true.

       nullify(current_buf)

    end if

    return

50  call error_handler("MolMech: input CORE_SHELL_LIST is wrong")
100 call error_handler("MolMech: input TWO_BODY_LIST is wrong")
200 call error_handler("MolMech: input THREE_BODY_LIST is wrong")
300 call error_handler("MolMech: input FOUR_BODY_LIST is wrong")
400 call error_handler("MolMech: input STR_STR_LIST is wrong")
500 call error_handler("MolMech: input BND_BND_LIST is wrong")
600 call error_handler("MolMech: input STR_BND_LIST is wrong")
700 call error_handler("MolMech: input STR_TRS_LIST is wrong")

  end function read_Nb_lists
  !******************************************************************

  !******************************************************************
  subroutine resorting_Nb_lists()

    integer(kind=i4_kind) :: i,j,k,l,length1,length2

    i_n_cs: do i=1,n_cs
       length1=size(core_shell(i)%list,1)
       length2=size(core_shell(i)%list,2)
       j_length20: do j=1,length2
          l_length10: do l=1,length1
             k_n_sp0: do k=1,n_species
                if(core_shell(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   core_shell(i)%list(l,j)=k
                   exit k_n_sp0
                end if
             end do k_n_sp0
          end do l_length10
       end do j_length20
    end do i_n_cs

    i_n_2b: do i=1,n_2b
       length1=size(two_body(i)%list,1)
       length2=size(two_body(i)%list,2)
       j_length2: do j=1,length2
          l_length1: do l=1,length1
             k_n_sp: do k=1,n_species
                if(two_body(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   two_body(i)%list(l,j)=k
                   exit k_n_sp
                end if
             end do k_n_sp
          end do l_length1
       end do j_length2
    end do i_n_2b

    i_n_3b: do i=1,n_3b
       length1=size(three_body(i)%list,1)
       length2=size(three_body(i)%list,2)
       j_length21: do j=1,length2
          l_length11: do l=1,length1
             k_n_sp1: do k=1,n_species
                if(three_body(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   three_body(i)%list(l,j)=k
                   exit k_n_sp1
                end if
             end do k_n_sp1
          end do l_length11
       end do j_length21
    end do i_n_3b

    i_n_4b: do i=1,n_4b
       length1=size(four_body(i)%list,1)
       length2=size(four_body(i)%list,2)
       j_length22: do j=1,length2
          l_length12: do l=1,length1
             k_n_sp2: do k=1,n_species
                if(four_body(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   four_body(i)%list(l,j)=k
                   exit k_n_sp2
                end if
             end do k_n_sp2
          end do l_length12
       end do j_length22
    end do i_n_4b

    i_n_ss: do i=1,n_ss
       length1=size(str_str(i)%list,1)
       length2=size(str_str(i)%list,2)
       j_length23: do j=1,length2
          l_length13: do l=1,length1
             k_n_sp3: do k=1,n_species
                if(str_str(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   str_str(i)%list(l,j)=k
                   exit k_n_sp3
                end if
             end do k_n_sp3
          end do l_length13
       end do j_length23
    end do i_n_ss

    i_n_bb: do i=1,n_bb
       length1=size(bnd_bnd(i)%list,1)
       length2=size(bnd_bnd(i)%list,2)
       j_length24: do j=1,length2
          l_length14: do l=1,length1
             k_n_sp4: do k=1,n_species
                if(bnd_bnd(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   bnd_bnd(i)%list(l,j)=k
                   exit k_n_sp4
                end if
             end do k_n_sp4
          end do l_length14
       end do j_length24
    end do i_n_bb

    i_n_sb: do i=1,n_sb
       length1=size(str_bnd(i)%list,1)
       length2=size(str_bnd(i)%list,2)
       j_length25: do j=1,length2
          l_length15: do l=1,length1
             k_n_sp5: do k=1,n_species
                if(str_bnd(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   str_bnd(i)%list(l,j)=k
                   exit k_n_sp5
                end if
             end do k_n_sp5
          end do l_length15
       end do j_length25
    end do i_n_sb

    i_n_st: do i=1,n_st
       length1=size(str_trs(i)%list,1)
       length2=size(str_trs(i)%list,2)
       j_length26: do j=1,length2
          l_length16: do l=1,length1
             k_n_sp6: do k=1,n_species
                if(str_trs(i)%list(l,j) == atoms_cart(k)%initial_number) then
                   str_trs(i)%list(l,j)=k
                   exit k_n_sp6
                end if
             end do k_n_sp6
          end do l_length16
       end do j_length26
    end do i_n_st

  end subroutine resorting_Nb_lists
  !******************************************************************

  !******************************************************************
  subroutine autobuilding_nb_lists()
    !automatic generation of N-body lists of bonded interactions
    !basing on atom connectivity calculated with help of covalent
    !radii of atoms (species) and core-shell cutoff
    !NOT COMPLETED (2-, 3- and c-s interactions, only)!!!!!
    
    integer(kind=i4_kind) :: i,j,k,l,m,n,kn,status,count,count1,nl
    integer(kind=i4_kind) :: type_buf,ty(10),in(10),i1,i2,i3,i4
    real(kind=r8_kind) :: rr_cov,r(3),dr,rr1,rr2,r1(3),r2(3),dr1,dr2,rc1,rc2
    real(kind=r8_kind) :: r3(3),r4(3),rc3,rc4
    real(r8_kind), allocatable :: buffer(:,:) 

    type buf_data
       integer(kind=i4_kind) :: ind(10)
       type (buf_data), pointer :: next_buf
    end type buf_data
    type (buf_data), target :: first_buf
    type (buf_data), pointer :: current_buf, tmp_buf,del

    !core-shell lists
    if(n_cs > 0) then
!!$print*,"CORE_SHELL_LIST"
       allocate(core_shell(n_cs),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed CORE_SHELL allocation")
       nl=0
       kn=n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+n_st+n_2b_n
       l1a: do i=kn+1,kn+n_cs
          current_buf=>first_buf
          nullify(current_buf%next_buf)

          count=0
          ty(1)=poten(i)%sp_type(1)
          ty(2)=poten(i)%sp_type(2)
          l3a: do j=1,n_species
             if(atoms_cart(j)%type /= ty(1) .and. atoms_cart(j)%type /= ty(2)) cycle l3a
             if(atoms_cart(j)%type == ty(1)) then 
                in(1)=1; in(2)=2
             end if
             if(atoms_cart(j)%type == ty(2)) then
                in(1)=2; in(2)=1
             end if
           
             l4a: do k=j+1,n_species
                if(atoms_cart(k)%type /= ty(in(2))) cycle l4a
                r=atoms_cart(j)%r-atoms_cart(k)%r
                if(lattice_calc) then 
                   call image(r)
                else if(slab_calc) then
                   call image_slab(r)
                end if
                dr=sqrt(dot_product(r,r))
                if(dr <= rcut_cs) then
                   allocate(tmp_buf)
                   tmp_buf%ind(1)=j; tmp_buf%ind(2)=k
                   nullify(tmp_buf%next_buf)
                   current_buf%next_buf=>tmp_buf
                   current_buf=>tmp_buf

                   count=count+1
                end if
             end do l4a
          end do l3a

          k=i-kn
          nl=nl+count
          allocate(core_shell(k)%list(2,count),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed CORE_SHELL%LIST allocation")
          current_buf=>first_buf
          do j=1,count
             del=>current_buf
             current_buf=>current_buf%next_buf
             core_shell(k)%list(1,j)=current_buf%ind(1)
             core_shell(k)%list(2,j)=current_buf%ind(2)
             if(j > 1) deallocate(del)
!!$print*,k,core_shell(k)%list(1,j),core_shell(k)%list(2,j)
          end do
       end do l1a

       if(to_epe) then
          allocate(cs_pair(2,nl),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed CS_PAIR allocation 1")
          k=0
          do i=1,n_cs
             l=size(core_shell(i)%list,2)
             do j=1,l
                k=k+1
                cs_pair(1,k)=core_shell(i)%list(1,j)
                cs_pair(2,k)=core_shell(i)%list(2,j)
             end do
          end do
       end if

    end if

    !2-body lists
!!$print*,n_2b,'n_2b'
    if(n_2b > 0) then
       allocate(two_body(n_2b),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed TWO_BODY allocation")
       l1: do i=1,n_2b
          current_buf=>first_buf
          nullify(current_buf%next_buf)

          count=0
          ty(1)=poten(i)%sp_type(1)
          ty(2)=poten(i)%sp_type(2)
 
          l3: do j=1,n_species
             if(atoms_cart(j)%type /= ty(1) .and. atoms_cart(j)%type /= ty(2)) cycle l3
             if(atoms_cart(j)%type == ty(1)) then 
                in(1)=1; in(2)=2
             end if
             if(atoms_cart(j)%type == ty(2)) then
                in(1)=2; in(2)=1
             end if
             i1=atoms_cart(j)%type
             if(ty(in(2)) <= n_species_types) then ! define 2-body list between regular species
                l4: do k=j+1,n_species
                   if(atoms_cart(k)%type /= ty(in(2))) cycle l4
                   i2=atoms_cart(k)%type
                   rr_cov=atoms(i1)%r_coval+atoms(i2)%r_coval
                   r=atoms_cart(j)%r-atoms_cart(k)%r
                   if(lattice_calc) then 
                      call image(r)
                   else if(slab_calc) then
                      call image_slab(r)
                   end if
                   dr=sqrt(dot_product(r,r))
                   if(dr <= rr_cov*rc_inc) then
                      if(poten(i)%id == 16) then         !for bck_short potential
                         ty(3)=poten(i)%sp_type(3)
                         if(ty(3) == 0) then
                            allocate(tmp_buf)
                            tmp_buf%ind(1)=j; tmp_buf%ind(2)=k
                            nullify(tmp_buf%next_buf)
                            current_buf%next_buf=>tmp_buf
                            current_buf=>tmp_buf
                            
                            count=count+1
                         else
                            do l=1,n_species
                               if(l==j .or. l==k) cycle
                               i3=atoms_cart(l)%type
                               if(i3 /= ty(3)) cycle
                               rr1=atoms(i1)%r_coval+atoms(i3)%r_coval
                               rr2=atoms(i2)%r_coval+atoms(i3)%r_coval
                               r1=atoms_cart(j)%r-atoms_cart(l)%r
                               if(lattice_calc) then 
                                  call image(r1)
                               else if(slab_calc) then
                                  call image_slab(r1)
                               end if
                               r2=atoms_cart(k)%r-atoms_cart(l)%r
                               if(lattice_calc) then 
                                  call image(r2)
                               else if(slab_calc) then
                                  call image_slab(r2)
                               end if
                               dr1=sqrt(dot_product(r1,r1))
                               dr2=sqrt(dot_product(r2,r2))
                               if(dr1 <= rr1 .and. dr2 <= rr2) then
                                  allocate(tmp_buf)
                                  tmp_buf%ind(1)=j; tmp_buf%ind(2)=k
                                  nullify(tmp_buf%next_buf)
                                  current_buf%next_buf=>tmp_buf
                                  current_buf=>tmp_buf
                                  
                                  count=count+1
                                  exit
                               end if
                            end do
                         end if
                      else                    ! usual potentials
                         allocate(tmp_buf)
                         tmp_buf%ind(1)=j; tmp_buf%ind(2)=k
                         nullify(tmp_buf%next_buf)
                         current_buf%next_buf=>tmp_buf
                         current_buf=>tmp_buf
                         
                         count=count+1
                      end if
                   end if
                end do l4
             else !define 2-body list between regular species and external PC
                l5: do k=1,n_ext_pc
                   if(pc_cart(k)%type /= ty(in(2))) cycle l5
                   i2=pc_cart(k)%type-n_species_types
                   rr_cov=atoms(i1)%r_coval+p_charges(i2)%r_coval
                   r=atoms_cart(j)%r-pc_cart(k)%r
                   dr=sqrt(dot_product(r,r))
                   if(dr <= rr_cov*rc_inc) then
                      allocate(tmp_buf)
                      tmp_buf%ind(1)=j; tmp_buf%ind(2)=k+n_species
                      nullify(tmp_buf%next_buf)
                      current_buf%next_buf=>tmp_buf
                      current_buf=>tmp_buf
                      
                      count=count+1
                   end if
                end do l5
             end if
          end do l3
          allocate(two_body(i)%list(2,count),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed TWO_BODY%LIST allocation")
          current_buf=>first_buf
          do j=1,count
             del=>current_buf
             current_buf=>current_buf%next_buf
             two_body(i)%list(1,j)=current_buf%ind(1)
             two_body(i)%list(2,j)=current_buf%ind(2)
             if(j > 1) deallocate(del)
!!$print*,i,two_body(i)%list(1,j),two_body(i)%list(2,j)
          end do
       end do l1
    end if

    !3-body lists
!!$print*,n_3b,'n_3b'
    if(n_3b > 0) then
       allocate(three_body(n_3b),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed THREE_BODY allocation")
       ll1: do i=n_2b+1,n_2b+n_3b
          current_buf=>first_buf
          nullify(current_buf%next_buf)

          count=0
          ty(1)=poten(i)%sp_type(1)
          ty(2)=poten(i)%sp_type(2)
          ty(3)=poten(i)%sp_type(3)

          ll3: do j=1,n_species+n_ext_pc
             if(j <= n_species) then 
                type_buf=atoms_cart(j)%type
                r1=atoms_cart(j)%r
             end if
             if(j > n_species) then 
                type_buf=pc_cart(j-n_species)%type
                r1=pc_cart(j-n_species)%r
             end if
             if(type_buf /= ty(2)) cycle ll3
             in(2)=2
             
             i1=type_buf
             if(i1 <= n_species_types) rc1=atoms(i1)%r_coval
             if(i1 > n_species_types) rc1=p_charges(i1-n_species_types)%r_coval
             ll4: do k=1,n_species+n_ext_pc
                if(k <= n_species) then 
                   type_buf=atoms_cart(k)%type
                   r2=atoms_cart(k)%r
                end if
                if(k > n_species) then 
                   type_buf=pc_cart(k-n_species)%type
                   r2=pc_cart(k-n_species)%r
                end if
                if(k==j) cycle ll4
                if(type_buf /= ty(1) .and. type_buf /= ty(3)) cycle ll4
                if(type_buf == ty(1)) then 
                   in(1)=1; in(3)=3
                end if
                if(type_buf == ty(3)) then 
                   in(1)=3; in(3)=1
                end if
                i2=type_buf
                if(i2 <= n_species_types) rc2=atoms(i2)%r_coval
                if(i2 > n_species_types) rc2=p_charges(i2-n_species_types)%r_coval
                rr_cov=rc1+rc2
                r=r1-r2
                if(lattice_calc) then 
                   call image(r)
                else if(slab_calc) then
                   call image_slab(r)
                end if
                dr=sqrt(dot_product(r,r))
                if(dr <= rr_cov*rc_inc) then
                   ll2: do l=1,n_species+n_ext_pc
                      if(l==j .or. l==k) cycle ll2
                      if(l < k) cycle ll2
                      if(l <= n_species) then
                         type_buf=atoms_cart(l)%type
                         r2=atoms_cart(l)%r
                      end if
                      if(l > n_species) then
                         type_buf=pc_cart(l-n_species)%type
                         r2=pc_cart(l-n_species)%r
                      end if
                      if(type_buf /= ty(in(3))) cycle ll2
                      i3=type_buf
                      if(i3 <= n_species_types) rc2=atoms(i3)%r_coval
                      if(i3 > n_species_types) rc2=p_charges(i3-n_species_types)%r_coval
                      rr_cov=rc1+rc2
                      r=r1-r2
                      if(lattice_calc) then 
                         call image(r)
                      else if(slab_calc) then
                         call image_slab(r)
                      end if
                      dr=sqrt(dot_product(r,r))
                      if(dr <= rr_cov*rc_inc) then
                         allocate(tmp_buf)
                         tmp_buf%ind(1)=k; tmp_buf%ind(2)=j; tmp_buf%ind(3)=l
                         nullify(tmp_buf%next_buf)
                         current_buf%next_buf=>tmp_buf
                         current_buf=>tmp_buf
                         
                         count=count+1
                      end if
                   end do ll2
                end if
             end do ll4
          end do ll3

          k=i-n_2b
          allocate(three_body(k)%list(3,count),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed THREE_BODY%LIST allocation")
          current_buf=>first_buf
          do j=1,count
             del=>current_buf
             current_buf=>current_buf%next_buf
             three_body(k)%list(1,j)=current_buf%ind(1)
             three_body(k)%list(2,j)=current_buf%ind(2)
             three_body(k)%list(3,j)=current_buf%ind(3)
             if(j > 1) deallocate(del)
!!$print*,k,three_body(k)%list(1,j),three_body(k)%list(2,j),three_body(k)%list(3,j)
          end do
       end do ll1
    end if

    !4-body lists
!!$print*,n_4b,'n_4b'
    if(n_4b > 0) then
       allocate(four_body(n_4b),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed FOUR_BODY allocation")
       kn=n_2b+n_3b
       lll1: do i=kn+1,kn+n_4b
          current_buf=>first_buf
          nullify(current_buf%next_buf)
          
          count=0
          ty(1)=poten(i)%sp_type(1)
          ty(2)=poten(i)%sp_type(2)
          ty(3)=poten(i)%sp_type(3)
          ty(4)=poten(i)%sp_type(4)
          
          lll2: do j=1,n_species
             type_buf=atoms_cart(j)%type
             if(type_buf /= ty(2) .and. type_buf /= ty(3)) cycle lll2
             if(type_buf == ty(2)) then
                in(1)=1; in(2)=2; in(3)=3; in(4)=4
             else
                cycle lll2
             end if
             r2=atoms_cart(j)%r
             i2=type_buf
             rc2=atoms(i2)%r_coval
             lll3: do k=1,n_species
                if(k==j) cycle lll3
                type_buf=atoms_cart(k)%type
                if(type_buf /= ty(in(3)))cycle lll3
                i3=type_buf
                rc3=atoms(i3)%r_coval
                r3=atoms_cart(k)%r
                r=r2-r3
                if(lattice_calc) then 
                   call image(r)
                else if(slab_calc) then
                   call image_slab(r)
                end if
                rr_cov=rc2+rc3
                dr=sqrt(dot_product(r,r))
                if(dr <= rr_cov*rc_inc) then
                   lll4: do l=1,n_species
                      if(l==j .or. l==k) cycle lll4
                      type_buf=atoms_cart(l)%type
                      if(type_buf /= ty(in(1))) cycle lll4
                      i1=type_buf
                      rc1=atoms(i1)%r_coval
                      r1=atoms_cart(l)%r
                      r=r2-r1
                      if(lattice_calc) then 
                         call image(r)
                      else if(slab_calc) then
                         call image_slab(r)
                      end if
                      rr_cov=rc2+rc1
                      dr=sqrt(dot_product(r,r))                   
                      if(dr <= rr_cov*rc_inc) then
                         lll5: do m=1,n_species
                            if(m==j .or. m==k .or. m==l) cycle lll5
                            type_buf=atoms_cart(m)%type
                            if(type_buf /= ty(in(4))) cycle lll5
                            i4=type_buf
                            rc4=atoms(i4)%r_coval
                            r4=atoms_cart(m)%r
                            r=r3-r4
                            if(lattice_calc) then 
                               call image(r)
                            else if(slab_calc) then
                               call image_slab(r)
                            end if
                            rr_cov=rc4+rc3
                            dr=sqrt(dot_product(r,r))
                            if(dr <= rr_cov*rc_inc) then
                               allocate(tmp_buf)
                               tmp_buf%ind(1)=l; tmp_buf%ind(2)=j; tmp_buf%ind(3)=k; tmp_buf%ind(4)=m
                               nullify(tmp_buf%next_buf)
                               current_buf%next_buf=>tmp_buf
                               current_buf=>tmp_buf
                               
                               count=count+1
                            end if
                         end do lll5
                      end if
                   end do lll4
                end if
             end do lll3
          end do lll2

          n=i-kn
          allocate(buffer(4,count),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed FOUR_BODY BUFFER allocation")
          current_buf=>first_buf
          count1=0
          j_count: do j=1,count
             del=>current_buf
             current_buf=>current_buf%next_buf
             do k=1,count1
                if(current_buf%ind(1) == buffer(4,k) .and. &
                    current_buf%ind(2) == buffer(3,k) .and. &
                    current_buf%ind(3) == buffer(2,k) .and. &
                    current_buf%ind(4) == buffer(1,k)) cycle j_count
             end do
             buffer(1,count1+1)=current_buf%ind(1)
             buffer(2,count1+1)=current_buf%ind(2)
             buffer(3,count1+1)=current_buf%ind(3)
             buffer(4,count1+1)=current_buf%ind(4)
             count1=count1+1
             if(j > 1) deallocate(del)
          end do j_count
          
          allocate(four_body(n)%list(4,count1),stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed FOUR_BODY%LIST allocation")
          four_body(n)%list(1:4,1:count1)=buffer(1:4,1:count1)

          deallocate(buffer,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed FOUR_BODY BUFFER deallocation")
          
!!$do j=1,count1
!!$print*,n,four_body(n)%list(1,j),four_body(n)%list(2,j),four_body(n)%list(3,j),four_body(n)%list(4,j)
!!$end do
       end do lll1
    end if

  end subroutine autobuilding_nb_lists
  !******************************************************************

  !******************************************************************
  subroutine build_2b_nonbonded_lists(generate_exclude_list,only_exclude_list)

    logical,intent(in) :: generate_exclude_list,only_exclude_list
    type buf_data
       integer(kind=i4_kind) :: ind(3)
       logical :: yes
       type (buf_data), pointer :: next_buf
    end type buf_data
    type (buf_data), target :: first_buf
    type (buf_data), pointer :: current_buf, tmp_buf, del

    integer(kind=i4_kind), allocatable :: buffer(:)
    integer(kind=i4_kind) :: i,j,k,l,status,length,count0,count1,count2,count3
    integer(kind=i4_kind) :: ia1,ia2,ia3,ia4,im1,im2,im3,jj,kk,ll,kn,ims,ift,n_t
    integer(kind=i4_kind) :: ty1,ty2
    logical :: yes,yess,yesss
    real(kind=r8_kind) :: dr,r(3)

    if(.not.coulomb .and. n_2b_n==0) return
 
    !definition atoms for full direct Coulomb summation and for 
    !building and rebuilding van der Waals 2-body lists
    if(.not.allocated(N_total)) then
       allocate(N_total(n_species),stat=status)
       if(status /= 0) call error_handler("MolMech: failed N_total allocation")
    end if

    if(mod(n_species,2)==0) then
       do i=1,n_species
          if(i<=n_species/2 .and. mod(i,2)==1) N_total(i)=n_species/2
          if(i<=n_species/2 .and. mod(i,2)==0) N_total(i)=n_species/2-1
          if(i >n_species/2 .and. mod(i-n_species/2,2)==1) N_total(i)=n_species/2-1
          if(i >n_species/2 .and. mod(i-n_species/2,2)==0) N_total(i)=n_species/2
       end do
    else if(mod(n_species-1,2)==0) then
       N_total=(n_species-1)/2
    end if

    if(generate_exclude_list) then 
       call exluded_list()
       if(lattice_calc .or.slab_calc) then 
          if(coulomb) call exluded_list_c()
       end if
       if (only_exclude_list) return
    end if
    
    if(n_2b_n/=0 .and. van_der_waals) call vdw_list()

    if(coulomb) then 
       if (lattice_calc) then 
          call ewald_list()
       else if(slab_calc) then
          call ewald2d_list()
       end if
    end if

  contains
    !-----------------------------------------------------------------
    subroutine exluded_list()
      !creating excluded lists for each species to avoid Coulomb and
      !van der Waals interactions between bonded (via bond or angle)
      !atoms

      allocate(bonded_spcs(n_species),stat=status)
      if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS allocation")
!!$print*,"EXCLUDE LIST"
      do i=1,n_species
         current_buf=>first_buf
         nullify(current_buf%next_buf)

         count0=0

         do j=1,n_cs
            length=size(core_shell(j)%list,2)

            do k=1,length
               im1=0
               ia1=core_shell(j)%list(1,k)
               ia2=core_shell(j)%list(2,k)
               if(i == ia1) im1=ia2 
               if(i == ia2) im1=ia1
               if(im1 /= 0) then
                  allocate(tmp_buf,stat=status)
                  if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl)")
                  tmp_buf%ind(1)=im1
                  nullify(tmp_buf%next_buf)
                  current_buf%next_buf=>tmp_buf
                  current_buf=>tmp_buf

                  count0=count0+1
               end if
            end do
         end do

         if(exclude_1_2) then
            do j=1,n_2b
               length=size(two_body(j)%list,2)

               do k=1,length
                  im1=0
                  ia1=two_body(j)%list(1,k)
                  ia2=two_body(j)%list(2,k)
                  if(i == ia1) im1=ia2 
                  if(i == ia2) im1=ia1
                  if(im1 /= 0) then
                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl1)")
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1
                  end if
               end do
            end do
         end if

         if(exclude_1_3) then
            do j=1,n_3b
               length=size(three_body(j)%list,2)

               do k=1,length
                  im1=0
                  im2=0
                  ia1=three_body(j)%list(1,k)
                  ia2=three_body(j)%list(2,k)
                  ia3=three_body(j)%list(3,k)
                  if(i == ia1) then
                     im1=ia2; im2=ia3 
                  else if(i == ia2) then
                     im1=ia1; im2=ia3
                  else if(i == ia3) then
                     im1=ia1; im2=ia2
                  end if
                  if(im1 /= 0 .and. im2 /= 0) then
                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl2)")
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl3)")
                     tmp_buf%ind(1)=im2
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1
                  end if
               end do
            end do

            do j=1,n_4b
               length=size(four_body(j)%list,2)
               
               do k=1,length
                  im1=0; im2=0; im3=0 
                  ia1=four_body(j)%list(1,k)
                  ia2=four_body(j)%list(2,k)
                  ia3=four_body(j)%list(3,k)
                  ia4=four_body(j)%list(4,k)
                  if(i==ia1) then
                     im1=ia2; im2=ia3
                  elseif(i==ia4) then
                     im1=ia2; im2=ia3
                  elseif(i==ia2) then
                     im1=ia1; im2=ia3; im3=ia4
                  elseif(i==ia3) then
                     im1=ia1; im2=ia2; im3=ia4
                  end if
                  if(im1 /= 0 .and. im2 /= 0) then
                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl4)")
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl5)")
                     tmp_buf%ind(1)=im2
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     if(im3 /= 0) then
                        allocate(tmp_buf,stat=status)
                        if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl6)")
                        tmp_buf%ind(1)=im3
                        nullify(tmp_buf%next_buf)
                        current_buf%next_buf=>tmp_buf
                        current_buf=>tmp_buf

                        count0=count0+1
                     end if
                  end if
               end do
            end do
         end if

         if(exclude_1_4) then
            do j=1,n_4b
               length=size(four_body(j)%list,2)
               
               do k=1,length
                  im1=0; im2=0; im3=0
                  ia1=four_body(j)%list(1,k)
                  ia2=four_body(j)%list(2,k)
                  ia3=four_body(j)%list(3,k)
                  ia4=four_body(j)%list(4,k)
                  if(i==ia1) then
                     im1=ia2; im2=ia3; im3=ia4
                  elseif(i==ia2) then
                     im1=ia1; im2=ia3; im3=ia4
                  elseif(i==ia3) then
                     im1=ia1; im2=ia2; im3=ia4
                  elseif(i==ia4) then
                     im1=ia1; im2=ia2; im3=ia3
                  end if
                  if(im1 /= 0 .and. im2 /= 0 .and. im3 /= 0) then
                     allocate(tmp_buf)
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf)
                     tmp_buf%ind(1)=im2
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf)
                     tmp_buf%ind(1)=im3
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1
                  end if
               end do
            end do
         end if
       
         allocate(buffer(count0),stat=status)
         if(status /= 0) call error_handler("MolMech: failed BUFFER allocation")
         current_buf=>first_buf
         count1=0
         j_count: do j=1,count0
            current_buf=>current_buf%next_buf
            do k=1,count1
               if(current_buf%ind(1) == buffer(k)) cycle j_count
            end do
            buffer(count1+1)=current_buf%ind(1)
            count1=count1+1
         end do j_count

         allocate(bonded_spcs(i)%list(count1,1), stat=status)
         if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS%LIST allocation")
         bonded_spcs(i)%list(1:count1,1)=buffer(1:count1)

         deallocate(buffer,stat=status)
         if(status /= 0) call error_handler("MolMech: failed BUFFER deallocation")

         do j=2,count1
            jj=bonded_spcs(i)%list(j,1)
            do k=1,j-1
               kk=bonded_spcs(i)%list(k,1)
               if(jj >= kk ) cycle
               ll=jj
               do l=j-1,k,-1
                  bonded_spcs(i)%list(l+1,1)=bonded_spcs(i)%list(l,1)
               end do
               bonded_spcs(i)%list(k,1)=ll
               exit
            end do
         end do

         nullify(current_buf)

!!$print*,'!',i,bonded_spcs(i)%list(:,1),size(bonded_spcs(i)%list,1)
      end do

    end subroutine exluded_list
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    subroutine exluded_list_c()
      !creating excluded lists for each species to avoid Coulomb (Ewald)
      !interaction between bonded (via bond or angle)
      !atoms

      allocate(bonded_spcs_c(n_species),stat=status)
      if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS_C allocation")
!!$print*,"EXCLUDE LIST"
      do i=1,n_species
         current_buf=>first_buf
         nullify(current_buf%next_buf)

         count0=0

         do j=1,n_cs
            length=size(core_shell(j)%list,2)

            do k=1,length
               im1=0
               ia1=core_shell(j)%list(1,k)
               ia2=core_shell(j)%list(2,k)
               if(i == ia1) im1=ia2
               if(i == ia2) im1=ia1
               if(im1 /= 0) then
                  allocate(tmp_buf,stat=status)
                  if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl_c)")
                  tmp_buf%ind(1)=im1
                  nullify(tmp_buf%next_buf)
                  current_buf%next_buf=>tmp_buf
                  current_buf=>tmp_buf

                  count0=count0+1
               end if
            end do
         end do

         if(exclude_1_2 .and. .not. ff3_style) then
            do j=1,n_2b
               length=size(two_body(j)%list,2)

               do k=1,length
                  im1=0
                  ia1=two_body(j)%list(1,k)
                  ia2=two_body(j)%list(2,k)
                  if(i == ia1) im1=ia2
                  if(i == ia2) im1=ia1
                  if(im1 /= 0) then
                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl1_c)")
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1
                  end if
               end do
            end do
         end if

         if(exclude_1_3) then
            do j=1,n_3b
               length=size(three_body(j)%list,2)

               do k=1,length
                  im1=0
                  im2=0
                  ia1=three_body(j)%list(1,k)
                  ia2=three_body(j)%list(2,k)
                  ia3=three_body(j)%list(3,k)
                  if(i == ia1) then
                     im1=ia2; im2=ia3
                  else if(i == ia2) then
                     im1=ia1; im2=ia3
                  else if(i == ia3) then
                     im1=ia1; im2=ia2
                  end if
                  if(im1 /= 0 .and. im2 /= 0) then
                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl2_c)")
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl3_c)")
                     tmp_buf%ind(1)=im2
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1
                  end if
               end do
            end do

            do j=1,n_4b
               length=size(four_body(j)%list,2)

               do k=1,length
                  im1=0; im2=0; im3=0
                  ia1=four_body(j)%list(1,k)
                  ia2=four_body(j)%list(2,k)
                  ia3=four_body(j)%list(3,k)
                  ia4=four_body(j)%list(4,k)
                  if(i==ia1) then
                     im1=ia2; im2=ia3
                  elseif(i==ia4) then
                     im1=ia2; im2=ia3
                  elseif(i==ia2) then
                     im1=ia1; im2=ia3; im3=ia4
                  elseif(i==ia3) then
                     im1=ia1; im2=ia2; im3=ia4
                  end if
                  if(im1 /= 0 .and. im2 /= 0) then
                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl4_c)")
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf,stat=status)
                     if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl5_c)")
                     tmp_buf%ind(1)=im2
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     if(im3 /= 0) then
                        allocate(tmp_buf,stat=status)
                        if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(excl6)")
                        tmp_buf%ind(1)=im3
                        nullify(tmp_buf%next_buf)
                        current_buf%next_buf=>tmp_buf
                        current_buf=>tmp_buf

                        count0=count0+1
                     end if
                  end if
               end do
            end do
         end if

         if(exclude_1_4) then
            do j=1,n_4b
               length=size(four_body(j)%list,2)

               do k=1,length
                  im1=0; im2=0; im3=0
                  ia1=four_body(j)%list(1,k)
                  ia2=four_body(j)%list(2,k)
                  ia3=four_body(j)%list(3,k)
                  ia4=four_body(j)%list(4,k)
                  if(i==ia1) then
                     im1=ia2; im2=ia3; im3=ia4
                  elseif(i==ia2) then
                     im1=ia1; im2=ia3; im3=ia4
                  elseif(i==ia3) then
                     im1=ia1; im2=ia2; im3=ia4
                  elseif(i==ia4) then
                     im1=ia1; im2=ia2; im3=ia3
                  end if
                  if(im1 /= 0 .and. im2 /= 0 .and. im3 /= 0) then
                     allocate(tmp_buf)
                     tmp_buf%ind(1)=im1
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf)
                     tmp_buf%ind(1)=im2
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1

                     allocate(tmp_buf)
                     tmp_buf%ind(1)=im3
                     nullify(tmp_buf%next_buf)
                     current_buf%next_buf=>tmp_buf
                     current_buf=>tmp_buf

                     count0=count0+1
                  end if
               end do
            end do
         end if

         allocate(buffer(count0),stat=status)
         if(status /= 0) call error_handler("MolMech: failed BUFFER allocation c")
         current_buf=>first_buf
         count1=0
         j_count: do j=1,count0
            current_buf=>current_buf%next_buf
            do k=1,count1
               if(current_buf%ind(1) == buffer(k)) cycle j_count
            end do
            buffer(count1+1)=current_buf%ind(1)
            count1=count1+1
         end do j_count

         allocate(bonded_spcs_c(i)%list(count1,1), stat=status)
         if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS_C%LIST allocation")
         bonded_spcs_c(i)%list(1:count1,1)=buffer(1:count1)

         deallocate(buffer,stat=status)
         if(status /= 0) call error_handler("MolMech: failed BUFFER deallocation c")

         do j=2,count1
            jj=bonded_spcs_c(i)%list(j,1)
            do k=1,j-1
               kk=bonded_spcs_c(i)%list(k,1)
               if(jj >= kk ) cycle
               ll=jj
               do l=j-1,k,-1
                  bonded_spcs_c(i)%list(l+1,1)=bonded_spcs_c(i)%list(l,1)
               end do
               bonded_spcs_c(i)%list(k,1)=ll
               exit
            end do
         end do

         nullify(current_buf)

!!$print*,'!',i,'-',bonded_spcs_c(i)%list(:,1),size(bonded_spcs(i)%list,1)
      end do

    end subroutine exluded_list_c
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    subroutine vdw_list()

      real(r8_kind) :: rc2,r1(3)
      integer(i4_kind) :: type_buf
      
      rc2=rcut_vdw*rcut_vdw

!!$print*,"VAN DER WAALS LIST",n_2b_n
      !building lists for short-range non bonded (van der Waals) interactions
      if(.not. allocated(non_bonded)) then
         allocate(non_bonded(n_species),stat=status)
         if(status /= 0) call error_handler("MolMech: failed NON_BONDED allocation")
         do i=1,n_species
            nullify(non_bonded(i)%list); nullify(non_bonded(i)%first_image)
         end do
      end if

      if(.not.minimal_image) then
         if(.not. allocated(add_list)) then
            allocate(add_list(n_species),stat=status)
            if(status /= 0) call error_handler("MolMech: failed ADD_LIST allocation")
            do i=1,n_species
               n_t=N_total(i); if(.not.minimal_image) n_t=N_total(i)+1
               allocate(add_list(i)%list(n_images,n_t),stat=status)
               if(status /= 0) call error_handler("MolMech: failed ADD_LIST%LIST allocation")
            end do
         end if
         do i=1,n_species
            add_list(i)%list=0
         end do
      end if

      i_sp: do i=1,n_species
         length=size(bonded_spcs(i)%list,1)

         current_buf=>first_buf
         nullify(current_buf%next_buf)

         yes=.false.
         j_nb: do j=1,n_2b_n
            kn=n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+n_st+j
            ty1=poten(kn)%sp_type(1)
            ty2=poten(kn)%sp_type(2)
            if(atoms_cart(i)%type == ty1 .or. atoms_cart(i)%type == ty2) then
               yes=.true.
               exit j_nb
            end if
         end do j_nb
         if(.not.yes) cycle i_sp

         count0=0; count3=0
         ift=1; if(.not.minimal_image) ift=0
         j_n: do j=ift,N_total(i)+n_ext_pc
            if(j <= N_total(i)) then !regular centers
               jj=i+j
               if(jj > n_species) jj=jj-n_species
               type_buf=atoms_cart(jj)%type
               r1=atoms_cart(jj)%r
            elseif(j > N_total(i)) then !external centers
               jj=j-N_total(i)+n_species
               type_buf=pc_cart(jj-n_species)%type
               r1=pc_cart(jj-n_species)%r
            end if

            yesss=.true.
            if(minimal_image) then
               k_l: do k=1,length
                  if(jj == bonded_spcs(i)%list(k,1)) then
                     if(minimal_image) then 
                        cycle j_n
                     else
                        yesss=.false.
                     end if
                  end if
               end do k_l
            end if

            k_n: do k=1,n_2b_n
               kn=n_2b+n_3b+n_4b+n_ss+n_bb+n_sb+n_st+k
               ty1=poten(kn)%sp_type(1)
               ty2=poten(kn)%sp_type(2)
               if(atoms_cart(i)%type /= ty1 .and. atoms_cart(i)%type /= ty2) then 
                  cycle k_n
               else if(type_buf /= ty1 .and. type_buf /= ty2) then 
                  cycle k_n
               else if(atoms_cart(i)%type == ty1 .and. type_buf /= ty2) then
                  cycle k_n
               else if(type_buf /= ty1 .and. atoms_cart(i)%type == ty2) then
                  cycle k_n
               end if

               if(minimal_image) then
                  r=r1-atoms_cart(i)%r
                  if(lattice_calc) then 
                     call image(r)
                  else if(slab_calc) then
                     call image_slab(r)
                  end if
                  dr=dot_product(r,r)
                  yess=.true.
                  if(dr > rc2) cycle j_n
               else
                  count3=count3+1
                  r=r1-atoms_cart(i)%r
                  dr=dot_product(r,r)
                  yess=.true.
                  if(dr > rc2) yess=.false.
                  if(i == jj .or. .not.yesss) yess=.false.
                  count2=0
                  ims_n: do ims=1,n_images
                     r=im_coor(jj)%r(:,ims)-atoms_cart(i)%r
                     dr=dot_product(r,r)
                     if(dr > rc2+small1) cycle ims_n
                     count2=count2+1
                     add_list(i)%list(count2,count3)=ims
                  end do ims_n
                  if(.not.yess .and. count2==0) then
                     count3=count3-1
                     cycle j_n
                  end if
               end if

               allocate(tmp_buf,stat=status)
               if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(vdw)")
               tmp_buf%ind(1)=jj
               tmp_buf%ind(2)=kn
               tmp_buf%yes=yess
               nullify(tmp_buf%next_buf)
               current_buf%next_buf=>tmp_buf
               current_buf=>tmp_buf
               
               count0=count0+1
               
               cycle j_n
            end do k_n
         end do j_n

         if(associated(non_bonded(i)%list)) then
            deallocate(non_bonded(i)%list,non_bonded(i)%first_image,stat=status)
            if(status /= 0) call error_handler("MolMech: failed NON_BONDED(I)%LIST deallocation(1)")
            nullify(non_bonded(i)%list); nullify(non_bonded(i)%first_image)
         end if
         allocate(non_bonded(i)%list(2,count0),stat=status)
         if(status /= 0) call error_handler("MolMech: failed NON_BONDED%LIST allocation")
         allocate(non_bonded(i)%first_image(count0),stat=status)
         if(status /= 0) call error_handler("MolMech: failed NON_BONDED%FIRST_IMAGE allocation")

         current_buf=>first_buf
         do j=1,count0
            del=>current_buf
            current_buf=>current_buf%next_buf
            non_bonded(i)%list(1,j)=current_buf%ind(1)
            non_bonded(i)%list(2,j)=current_buf%ind(2)
            non_bonded(i)%first_image(j)=current_buf%yes
            if(j > 1) deallocate(del)
         end do

         nullify(current_buf)
!!$print*,'!',i,'|',size(non_bonded(i)%list,2),'-',non_bonded(i)%list(1,:)
      end do i_sp

    end subroutine vdw_list
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    subroutine ewald_list()

!!$print*,"EWALD LIST"
      !building lists for short-range non bonded (3D electrostatic) interactions
      if(.not. allocated(coulombic)) then
         allocate(coulombic(n_species),stat=status)
         if(status /= 0) call error_handler("MolMech: failed COULOMBIC allocation")
      end if

      if(.not.minimal_image) then
         if(.not. allocated(add_list_c)) then
            allocate(add_list_c(n_species),stat=status)
            if(status /= 0) call error_handler("MolMech: failed ADD_LIST_C allocation")
            do i=1,n_species
               n_t=N_total(i); if(.not.minimal_image) n_t=N_total(i)+1
               allocate(add_list_c(i)%list(n_images,n_t),stat=status)
               if(status /= 0) call error_handler("MolMech: failed ADD_LIST_C%LIST allocation")
            end do
         end if
         do i=1,n_species
            add_list_c(i)%list=0
         end do
      end if

      i_sp: do i=1,n_species
         length=size(bonded_spcs_c(i)%list,1)

         current_buf=>first_buf
         nullify(current_buf%next_buf)

         count0=0; count3=0
         ift=1; if(.not.minimal_image) ift=0
         j_n: do j=ift,N_total(i)
            jj=i+j
            if(jj > n_species) jj=jj-n_species

            yesss=.true.
!!$            k_l: do k=1,length
!!$               if(jj == bonded_spcs_c(i)%list(k,1)) then
!!$                  if(minimal_image) then 
!!$                     cycle j_n
!!$                  else
!!$                     yesss=.false.
!!$                  end if
!!$               end if
!!$            end do k_l

            if(minimal_image) then
               r=atoms_cart(jj)%r-atoms_cart(i)%r
               if(lattice_calc) call image(r)
               dr=sqrt(dot_product(r,r))
               yess=.true.
               if(dr > rcut_ew) cycle j_n
            else
               count3=count3+1
               r=atoms_cart(jj)%r-atoms_cart(i)%r
               dr=sqrt(dot_product(r,r))
               yess=.true.
               if(dr > rcut_ew) yess=.false.
               if(i == jj .or. .not.yesss) yess=.false.
               count2=0
               ims_n: do ims=1,n_images
                  r=im_coor(jj)%r(:,ims)-atoms_cart(i)%r
                  dr=sqrt(dot_product(r,r))
                  if(dr > rcut_ew+small) cycle ims_n
                  count2=count2+1
                  add_list_c(i)%list(count2,count3)=ims
               end do ims_n
               if(.not.yess .and. count2==0) then
                  count3=count3-1
                  cycle j_n
               end if
            end if
            
            allocate(tmp_buf,stat=status)
            if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(ew3list)")
            tmp_buf%ind(1)=jj
            tmp_buf%yes=yess
            nullify(tmp_buf%next_buf)
            current_buf%next_buf=>tmp_buf
            current_buf=>tmp_buf

            count0=count0+1

            cycle j_n
         end do j_n

         if(associated(coulombic(i)%list)) then
            deallocate(coulombic(i)%list,coulombic(i)%first_image,stat=status)
            if(status /= 0) call error_handler("MolMech: failed COULOMBI(I)%LIST deallocation(1)")
            nullify(coulombic(i)%list); nullify(coulombic(i)%first_image)
         end if
         allocate(coulombic(i)%list(count0,1),stat=status)
         if(status /= 0) call error_handler("MolMech: failed COULOMBIC%LIST allocation")
         allocate(coulombic(i)%first_image(count0),stat=status)
         if(status /= 0) call error_handler("MolMech: failed COULOMBIC%FIRST_IMAGE allocation")

         current_buf=>first_buf
         do j=1,count0
            del=>current_buf
            current_buf=>current_buf%next_buf
            coulombic(i)%list(j,1)=current_buf%ind(1)
            coulombic(i)%first_image(j)=current_buf%yes
            if(j > 1) deallocate(del)
         end do

         nullify(current_buf)
!!$print*,'!',i,coulombic(i)%list(1,:)
      end do i_sp

    end subroutine ewald_list
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    subroutine ewald2d_list()

!!$print*,"EWALD2D LIST"
      !building lists for short-range non bonded (2D electrostatic) interactions
      if(.not. allocated(coulombic)) then
         allocate(coulombic(n_species),stat=status)
         if(status /= 0) call error_handler("MolMech: failed COULOMBIC2D allocation")
      end if

      if(.not.minimal_image) then
         if(.not. allocated(add_list_c)) then
            allocate(add_list_c(n_species),stat=status)
            if(status /= 0) call error_handler("MolMech: failed ADD_LIST_C2D allocation")
            do i=1,n_species
               n_t=N_total(i); if(.not.minimal_image) n_t=N_total(i)+1
               allocate(add_list_c(i)%list(n_images,n_t),stat=status)
               if(status /= 0) call error_handler("MolMech: failed ADD_LIST_C%LIST2D allocation")
            end do
         end if
         do i=1,n_species
            add_list_c(i)%list=0
         end do
      end if

      i_sp: do i=1,n_species
         length=size(bonded_spcs_c(i)%list,1)

         current_buf=>first_buf
         nullify(current_buf%next_buf)

         count0=0; count3=0
         ift=1; if(.not.minimal_image) ift=0
         j_n: do j=ift,N_total(i)
            jj=i+j
            if(jj > n_species) jj=jj-n_species

            yesss=.true.
!!$            k_l: do k=1,length
!!$               if(jj == bonded_spcs_c(i)%list(k,1)) then
!!$                  if(minimal_image) then 
!!$                     cycle j_n
!!$                  else
!!$                     yesss=.false.
!!$                  end if
!!$               end if
!!$            end do k_l

            if(minimal_image) then
               r=atoms_cart(jj)%r-atoms_cart(i)%r
               if(slab_calc) call image_slab(r)
               dr=sqrt(dot_product(r,r))
               yess=.true.
               if(dr > rcut_ew) cycle j_n
            else
               count3=count3+1
               r=atoms_cart(jj)%r-atoms_cart(i)%r
               dr=sqrt(dot_product(r,r))
               yess=.true.
               if(dr > rcut_ew) yess=.false.
               if(i == jj .or. .not.yesss) yess=.false.
               count2=0
               ims_n: do ims=1,n_images
                  r=im_coor(jj)%r(:,ims)-atoms_cart(i)%r
                  dr=sqrt(dot_product(r,r))
                  if(dr > rcut_ew+small) cycle ims_n
                  count2=count2+1
                  add_list_c(i)%list(count2,count3)=ims
               end do ims_n
               if(.not.yess .and. count2==0) then
                  count3=count3-1
                  cycle j_n
               end if
            end if
            
            allocate(tmp_buf,stat=status)
            if(status /= 0) call error_handler("MolMech: failed tmp_data allocation(ew2list)")
            tmp_buf%ind(1)=jj
            tmp_buf%yes=yess
            nullify(tmp_buf%next_buf)
            current_buf%next_buf=>tmp_buf
            current_buf=>tmp_buf

            count0=count0+1

            cycle j_n
         end do j_n

         allocate(coulombic(i)%list(count0,1),stat=status)
         if(status /= 0) call error_handler("MolMech: failed COULOMBIC%LIST2D allocation")
         allocate(coulombic(i)%first_image(count0),stat=status)
         if(status /= 0) call error_handler("MolMech: failed COULOMBIC%FIRST_IMAGE2D allocation")

         current_buf=>first_buf
         do j=1,count0
            del=>current_buf
            current_buf=>current_buf%next_buf
            coulombic(i)%list(j,1)=current_buf%ind(1)
            coulombic(i)%first_image(j)=current_buf%yes
            if(j > 1) deallocate(del)
         end do

         nullify(current_buf)
!!$print*,'!',i,coulombic(i)%list(1,:)
      end do i_sp

    end subroutine ewald2d_list

  end subroutine build_2b_nonbonded_lists
  !******************************************************************

  !******************************************************************
  subroutine shutdown_interacting_lists()

    integer(i4_kind) :: status,i

    if(allocated(core_shell)) then
       do i=1,n_cs
          deallocate(core_shell(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed TWO_BODY(I)%LIST deallocation")
       end do
       deallocate(core_shell,stat=status)
       if(status /= 0) call error_handler("MolMech: failed TWO_BODY deallocation")
       if(allocated(cs_pair)) then
          deallocate(cs_pair,stat=status)
          if(status /= 0) call error_handler("MolMech: failed CS_PAIR deallocation")
       end if
    end if

    if(allocated(two_body)) then
       do i=1,n_2b
          deallocate(two_body(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed TWO_BODY(I)%LIST deallocation")
       end do
       deallocate(two_body,stat=status)
       if(status /= 0) call error_handler("MolMech: failed TWO_BODY deallocation")
    end if

    if(allocated(three_body)) then
       do i=1,n_3b
          deallocate(three_body(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed THREE_BODY(I)%LIST deallocation")
       end do
       deallocate(three_body,stat=status)
       if(status /= 0) call error_handler("MolMech: failed THREE_BODY deallocation")
    end if

    if(allocated(four_body)) then
       do i=1,n_4b
          deallocate(four_body(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed FOUR_BODY(I)%LIST deallocation")
       end do
       deallocate(four_body,stat=status)
       if(status /= 0) call error_handler("MolMech: failed FOUR_BODY deallocation")
    end if

    if(allocated(str_str)) then
       do i=1,n_ss
          deallocate(str_str(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed STR_STR(I)%LIST deallocation")
       end do
       deallocate(str_str,stat=status)
       if(status /= 0) call error_handler("MolMech: failed STR_STR deallocation")
    end if

    if(allocated(bnd_bnd)) then
       do i=1,n_bb
          deallocate(bnd_bnd(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed BND_BND(I)%LIST deallocation")
       end do
       deallocate(bnd_bnd,stat=status)
       if(status /= 0) call error_handler("MolMech: failed BND_BND deallocation")
    end if

    if(allocated(str_bnd)) then
       do i=1,n_sb
          deallocate(str_bnd(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed STR_BND(I)%LIST deallocation")
       end do
       deallocate(str_bnd,stat=status)
       if(status /= 0) call error_handler("MolMech: failed STR_BND deallocation")
    end if

    if(allocated(str_trs)) then
       do i=1,n_st
          deallocate(str_trs(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed STR_TRS(I)%LIST deallocation")
       end do
       deallocate(str_trs,stat=status)
       if(status /= 0) call error_handler("MolMech: failed STR_TRS deallocation")
    end if

    if(allocated(bonded_spcs)) then
       do i=1,n_species
          deallocate(bonded_spcs(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS(I)%LIST deallocation")
       end do
       deallocate(bonded_spcs,stat=status)
       if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS deallocation")
    end if

    if(allocated(bonded_spcs_c)) then
       do i=1,n_species
          deallocate(bonded_spcs_c(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS(I)%LIST deallocation")
       end do
       deallocate(bonded_spcs_c,stat=status)
       if(status /= 0) call error_handler("MolMech: failed BONDED_SPCS deallocation")
    end if

    if(allocated(non_bonded)) then
       do i=1,n_species
          if(associated(non_bonded(i)%list)) then
             deallocate(non_bonded(i)%list,non_bonded(i)%first_image,stat=status)
             if(status /= 0) call error_handler("MolMech: failed NON_BONDED(I)%LIST deallocation")
          end if
       end do
       deallocate(non_bonded,stat=status)
       if(status /= 0) call error_handler("MolMech: failed NON_BONDED deallocation")
    end if

    if(allocated(add_list)) then
       do i=1,n_species
          deallocate(add_list(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed ADD_LIST%LIST allocation")
       end do
       deallocate(add_list,stat=status)
       if(status /= 0) call error_handler("MolMech: failed ADD_LIST deallocation")
    end if

    if(allocated(coulombic)) then
       do i=1,n_species
          deallocate(coulombic(i)%list,coulombic(i)%first_image,stat=status)
          if(status /= 0) call error_handler("MolMech: failed COULOMBIC(I)%LIST deallocation")
       end do
       deallocate(coulombic,stat=status)
       if(status /= 0) call error_handler("MolMech: failed COULOMBIC deallocation")
    end if

    if(allocated(add_list_c)) then
       do i=1,n_species
          deallocate(add_list_c(i)%list,stat=status)
          if(status /= 0) call error_handler("MolMech: failed ADD_LIST_C%LIST allocation")
       end do
       deallocate(add_list_c,stat=status)
       if(status /= 0) call error_handler("MolMech: failed ADD_LIST_C deallocation")
    end if

    if(allocated(N_total)) then
       deallocate(N_total,stat=status)
       if(status /= 0) call error_handler("MolMech: failed N_TOTAL deallocation")
    end if

  end subroutine shutdown_interacting_lists
  !******************************************************************

  !******************************************************************
  subroutine send_receive_n_total()

    use molmech_msgtag_module
    use comm_module

    integer(i4_kind) :: info,status

    if(comm_i_am_master()) then
       call comm_init_send(comm_all_other_hosts,msgtag_mm_send_ntotal)
       
       call commpack(N_total(1),n_species,1,info)
       if( info /= 0) call error_handler &
            ("send_receive_n_total: n_total pack failed")

       call comm_send()
    else
       if(.not.allocated(N_total)) then 
          allocate(N_total(n_species),stat=status)
          if(status /= 0) call error_handler("MolMech: failed N_TOTAL allocation(1)")
       end if

       call communpack(N_total(1),n_species,1,info)
       if( info /= 0) call error_handler &
            ("send_receive_n_total: n_total unpack failed")
    end if

  end subroutine send_receive_n_total
  !******************************************************************

  !******************************************************************
end module n_body_lists_module






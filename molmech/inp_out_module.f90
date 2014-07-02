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
module inp_out_module

  !------------ Modules used -----------------------------------------
  use type_module
  use iounitadmin_module
  use filename_module

  implicit none
  private         ! by default, all names are private
  save

  !== Interrupt end of public interface of module ====================

  !------------ Declaration of public constants and variables -----
  integer(kind=i4_kind), public :: input_device, output_device

  !------------ public functions and subroutines ---------------------
  public get_input_device, get_output_device, close_input_device, &
       close_output_device, upcase, check_string, go_to_first_input_line, &
       name_without_cs, input_nm_error,repeated_definition,find_word,find_namelist, &
       get_file_device, close_file_device

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of private constants and variables ----
  integer(i4_kind) :: num
  !------------ Subroutines -----------------------------------------

contains
  !******************************************************************
  subroutine get_input_device(num1)

    integer(i4_kind), intent(in) :: num1
    logical :: do_preproc
    character(len=200) :: subfile
    character(len=250) :: message
    logical :: yes

    num=num1

    if(num==0) inquire (file=trim(inpfile('molmech.inp')), exist=yes)
    if(num==-1) inquire (file=trim(inpfile('molmech.inp')), exist=yes)
    if(num==2) inquire (file=trim(inpfile('molmech.inp.2')), exist=yes)
    if(num==3) inquire (file=trim(inpfile('molmech.inp.3')), exist=yes)
    if(yes) then
       do_preproc=preprocessing(input_device,subfile)
       if(.not.do_preproc) then
          message="Error occered during preprocessinf "//trim(subfile)//" input subfile"
          call error_handler(message)
       end if
    else
       call error_handler("MolMech: I cannot find the main input file (molmech.inp)")
    end if

  end subroutine  get_input_device
  !******************************************************************

  !******************************************************************
  function preprocessing(input_interm,subfile)

    logical :: preprocessing
    integer(kind=i4_kind), intent(out) :: input_interm
    character(len=200), intent(out) :: subfile

    integer(kind=i4_kind) :: i
    integer(kind=i4_kind) :: cur_file, file_id(100)
    integer(kind=i4_kind) :: ind
    character(len=200) :: input_line,input_line_save
    character(len=200) :: subfile_name(100)
    character(len=3) :: buf
    character(len=6) :: form
    logical :: yes

    input_interm=openget_iounit(trim(inpfile('molmech.inp.interm')), &
         form='formatted', status='unknown')

    cur_file=1
    if(num==0) subfile_name(1)=trim(inpfile('molmech.inp'))
    if(num==-1) subfile_name(1)=trim(inpfile('molmech.inp'))
    if(num==2) subfile_name(1)=trim(inpfile('molmech.inp.2'))
    if(num==3) subfile_name(1)=trim(inpfile('molmech.inp.3'))
    file_id(cur_file)=openget_iounit(subfile_name(1), form='formatted', status='old')

    main_cycle:do
       read(file_id(cur_file),'(a200)',end=100,err=200) input_line

       ind=index(input_line,'#')
       if (ind == 1) then
          input_line_save=input_line
          call upcase(input_line_save)
          if (check_string(input_line_save,"#INCLUDE ")) then
             cur_file=cur_file+1
             file_id(cur_file)=open_subfile()
             if (file_id(cur_file) == 0) then
                cur_file=cur_file-1
             end if
             cycle main_cycle
          else if(check_string(input_line_save,"#END")) then
             goto 100
          else
             cycle main_cycle
          end if
       end if
       if (ind == 0) ind=len_trim(input_line)+1
       if (ind == 1) cycle
       yes=.false.
       aa: do i=1,ind-1
          if(input_line(i:i) /=" ") then
             yes =.true.
             exit aa
          end if
       end do aa
       if(.not.yes) cycle
       write(buf,'(i3)') ind-1
       form='(a'//trim(adjustl(buf))//')'
       write(input_interm,form,err=200) trim(input_line)
       cycle main_cycle

100    call returnclose_iounit(file_id(cur_file))
       cur_file=cur_file-1
       if (cur_file == 0) exit main_cycle
    end do main_cycle

    preprocessing= .true.
    return

200 call returnclose_iounit(file_id(cur_file))
    preprocessing= .false.
    subfile=subfile_name(cur_file)
    return

  contains

    function open_subfile()

      integer(kind=i4_kind) :: open_subfile

      integer(kind=i4_kind) :: i
      character(len=200) :: line_buf
      logical :: yes

      i=9
      do
         if (input_line(i:i) == " ") then
            i=i+1
         else
            exit
         end if
      end do
      line_buf=input_line(i:200)
      if (check_string(line_buf,"/")) then
         subfile_name(cur_file)=trim(line_buf)
      else
         subfile_name(cur_file)=trim(inpfile(line_buf))
      end if

      inquire(file=subfile_name(cur_file),exist=yes)
      if (yes) then
         open_subfile=openget_iounit(subfile_name(cur_file), form='formatted', status='old')
      else
         open_subfile=0
      end if

    end function open_subfile
  end function preprocessing
  !******************************************************************

  !******************************************************************
  subroutine go_to_first_input_line()

    rewind input_device

  end subroutine go_to_first_input_line
  !******************************************************************

  !******************************************************************
  subroutine close_input_device()

    call returnclose_iounit(input_device,'delete')

  end subroutine close_input_device
  !******************************************************************

  !******************************************************************
  subroutine get_output_device(num1)

    integer(i4_kind), intent(in) :: num1

    if(num1 == 0) then
       output_device=openget_iounit(trim(outfile('molmech.out')), &
            form='formatted', status='unknown')
    end if
    if(num1 == -1) then
       output_device=openget_iounit(trim(outfile('molmech.out')), &
            form='formatted', status='unknown')
    end if
    if(num1 == 2) then
       output_device=openget_iounit(trim(outfile('molmech.out.2')), &
            form='formatted', status='unknown')
    end if
    if(num1 == 3) then
       output_device=openget_iounit(trim(outfile('molmech.out.3')), &
            form='formatted', status='unknown')
    end if

  end subroutine  get_output_device
  !******************************************************************

  !******************************************************************
  subroutine close_output_device()

    call returnclose_iounit(output_device)

  end subroutine close_output_device
  !******************************************************************

  !******************************************************************
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
  !******************************************************************

  !******************************************************************
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
  !******************************************************************

  !******************************************************************
  function find_word(word,i)

    character(len=*), intent(in) :: word
    integer(kind=i4_kind), intent(out) :: i

    logical :: find_word

    character(len=80) :: string, message
    character(len=4) :: number

    call go_to_first_input_line()
    i=0
    do
       i=i+1
       read(input_device,'(a80)', end=100, err=200) string
       call upcase(string)
       if (check_string(string,word)) then
          find_word= .true.
          exit
       end if
    end do
    return

100 find_word= .false.
    return
200 write(number,'(i4)') i
    message = trim(" MolMech: The error reading in line "//trim(number)// &
         " of intermediate input file")
    call error_handler(message)

  end function find_word
  !******************************************************************

  !******************************************************************
  function find_namelist(word,i)

    character(len=*), intent(in) :: word
    integer(kind=i4_kind), intent(out) :: i

    logical :: find_namelist

    character(len=80) :: string, message
    character(len=4) :: number

    i=0
    do
       i=i+1
       read(input_device,'(a80)', end=100, err=200) string
       call upcase(string)
       if (check_string(string,word)) then
          find_namelist= .true.
          exit
       end if
    end do
    backspace input_device
    return

100 find_namelist= .false.
    return
200 write(number,'(i4)') i
    message = trim(" MolMech: The error reading in line "//trim(number)// &
         " of intermediate input file")
    call error_handler(message)

  end function find_namelist
  !******************************************************************

  !******************************************************************
  subroutine name_without_cs(name_in,name_out)

    character(len=*), intent(in) :: name_in
    character(len=*), intent(out) :: name_out

    integer(kind=i4_kind) :: pos

    pos=index(name_in," c")
    if(pos /= 0_i4_kind) goto 1
    pos=index(name_in," C")
    if(pos /= 0_i4_kind) goto 1
    pos=index(name_in," s")
    if(pos /= 0_i4_kind) goto 1
    pos=index(name_in," S")
    if(pos /= 0_i4_kind) goto 1

    name_out=name_in
    return

1   name_out=name_in(1:pos-1)

  end subroutine name_without_cs
  !******************************************************************

  !******************************************************************
  subroutine input_nm_error(count,nm_name)

    integer(kind=i4_kind), intent(in) :: count
    character(len=*) :: nm_name

    character(len=4) :: number
    character(len=100) :: message

    if(count == 0) then
       number=" "
    else
       write(number,'(i4)') count
    end if
    message = trim(" MolMech: The error in "// &
         trim(number)//" namelist "//trim(nm_name))
    call error_handler(message)

  end subroutine input_nm_error
  !******************************************************************

  !******************************************************************
  subroutine repeated_definition(potential,name1,name2,name3,name4)

    character(len=*), intent(in) :: potential
    character(len=*), intent(in) :: name1,name2
    character(len=*), intent(in), optional :: name3,name4

    character(len=100) :: message

    if(.not.present(name3) .and. .not.present(name4)) then
       message = trim(" MolMech: Repeated definition of "// &
            trim(potential)//" interaction between "// &
            trim(name1)//" and "//trim(name2)//" atoms")
    else if(present(name3) .and. .not.present(name4)) then
       message = trim(" MolMech: Repeated definition of "// &
            trim(potential)//" interaction between "// &
            trim(name1)//", "//trim(name2)//" and"// &
            trim(name3)//" atoms")
    else if(present(name4)) then
       message = trim(" MolMech: Repeated definition of "// &
            trim(potential)//" interaction between "// &
            trim(name1)//", "//trim(name2)//", "// &
            trim(name3)//" and "//trim(name4)//" atoms")
    end if

    call error_handler(message)

  end subroutine repeated_definition
  !******************************************************************

  !******************************************************************
  subroutine get_file_device(file_device,file_name,inp_out)

    integer, intent(out) :: file_device
    character(len=*), intent(in) :: file_name
    character(len=3), intent(in) :: inp_out

    character(len=3) :: buf

    buf=inp_out

    call upcase(buf)

    if(buf == 'INP') then
       file_device=openget_iounit(trim(inpfile(file_name)), &
            form='formatted', status='unknown')
    elseif(buf == 'OUT') then
       file_device=openget_iounit(trim(outfile(file_name)), &
            form='formatted', status='unknown')
    else
       file_device=openget_iounit(trim(inpfile(file_name)), &
            form='formatted', status='unknown')
    end if

  end subroutine  get_file_device
  !******************************************************************

  !******************************************************************
  subroutine close_file_device(file_device)

    integer, intent(in) :: file_device

    call returnclose_iounit(file_device)

  end subroutine close_file_device
  !******************************************************************

  !******************************************************************

end module inp_out_module




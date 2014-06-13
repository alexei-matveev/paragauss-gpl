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
module qm_epe_interface_module
!-----------------------------------------------
#include <def.h>
  use type_module
  use datatype
  use filename_module, only: env => filename_env, filename_namelengthmax
  use iounitadmin_module
  use inp_out_module, only: upcase,check_string !from MM

  implicit none
  save        
  private

  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  type, public :: ff_param
     real(r8_kind) :: b
     real(r8_kind) :: r
     real(r8_kind) :: c
     real(r8_kind) :: d
     real(r8_kind) :: cut
  end type ff_param

  type(ff_param), allocatable, public :: interface_ff(:,:)

  integer(i4_kind), allocatable, public :: at_type(:)

  !------------ public functions and subroutines ------------------
   public read_interface_param, free_ff_arrays
   public find_block, read_block_line, check_comment, check_empty_line, wrong_line

  !================================================================
  ! End of public interface of module
  !================================================================

  integer(i4_kind), public :: line_number
  character(len=80), public :: current_block
  logical, public :: read_block,end_block,comment_line,empty_line
  character(len=10), parameter, public :: files(2)=(/"epeinput  ","forcefield"/)
  integer(i4_kind), parameter, public :: epeinp=1, fflib=2

  character(len=30), public :: ff_name
  integer(i4_kind), public :: ff_unit
  logical, public :: l_ff
  !------------ Declaration of private constants and variables ----

  integer(i4_kind) :: input_epe


  integer(i4_kind) :: n_type_atm

  !------------ Subroutines ---------------------------------------
contains
!====================================================================================

!====================================================================================
  function find_block(file_id,block_name,file_num)

    logical :: find_block
    integer(i4_kind), intent(in) :: file_id,file_num
    character(len=*), intent(in) :: block_name

    character(len=80) :: string

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

100 find_block= .false.
    return

200 call wrong_line(file_num)

  end function find_block
!====================================================================================

!====================================================================================
  subroutine read_block_line(file_id,line,end,comment,empty,file_num)

    integer(i4_kind), intent(in) :: file_id,file_num
    character(len=*), intent(out) :: line
    logical, intent(out) :: end,comment,empty

    character(len=80) :: string,message

    end=.false.; comment=.false.;empty=.false.
    line_number=line_number+1
    read(file_id,'(a80)', end=100, err=200) line
    string=line
    call upcase(string)
    if(check_empty_line(string)) empty=.true.
    if(check_comment(string)) comment=.true.
    if(check_string(string,"&END")) end=.true.
    
    return

100  message = trim("Read_"//files(file_num)//" : Unterminated block "//trim(current_block))
    call error_handler(message)

200 call wrong_line(file_num)

  end subroutine read_block_line
!====================================================================================

!====================================================================================
  function check_comment(string)

    logical :: check_comment
    character(len=*), intent(in) :: string

    check_comment=.false.
    if(index(string,"#") == 1) check_comment=.true.

  end function check_comment
!====================================================================================

!====================================================================================
  function check_empty_line(string)
    logical :: check_empty_line
    character(len=*), intent(in) :: string

    check_empty_line=.false.
    if(len(trim(string)) == 0) check_empty_line=.true.

  end function check_empty_line
!====================================================================================

!====================================================================================
  subroutine wrong_line(file_num)

    integer(i4_kind) :: file_num
    character(len=80) :: message
    character(len=4) :: number

    write(number,'(i4)') line_number
    message = trim("Read_"//files(file_num)//" : The error reading line "//trim(number))
    call error_handler(message)

  end subroutine wrong_line

!##############################################################################

!##############################################################################
  subroutine read_interface_param(n_atm,atnm,n_cluster)

    use slspar_module

    integer(i4_kind), intent(in) :: n_atm,n_cluster
    real(r8_kind), intent(in) :: atnm(:)

    character(len=filename_namelengthmax) :: libdir,input_dir
    character(len=80) :: buffer
    integer(i4_kind) :: i
    real(r8_kind) :: type_nm(100)

    if(.not.env("TTFSINPUTDIR",input_dir)) call error_handler( &
         "Read_interface_param: environment variable TTFSINPUTDIR is not set")

    inquire (file=trim(input_dir)//"/"//'/epe.input',exist=l_ff)
    if(l_ff) then 
       input_epe=openget_iounit(trim(input_dir)//'/epe.input', &
            form='formatted', status='old')
    else
       call error_handler("Read_interface_param: file "//trim(input_dir)//"/"//"epe.input"//" not found")
    end if

    read_block=find_block(input_epe,"&FORCE_FIELD",epeinp)
    if(.not.read_block) then
       call error_handler("Read_interface_param: No force field found within epe.input")
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

    if(.not.env("TTFSLIBS",libdir)) call error_handler( &
         "Read_interface_param: environment variable TTFSLIBS is not set")

    ff_name=adjustl(ff_name)
    inquire (file=trim(libdir)//"/"//trim(ff_name),exist=l_ff)
    if(l_ff) then 
       ff_unit=openget_iounit(trim(libdir)//"/"//trim(ff_name),  &
            form='formatted', status='unknown',action='read')
    else
       call error_handler("Read_interface_param: file "//trim(libdir)//"/"//trim(ff_name)//" not found")
    end if

    call get_type()

    call read_2body_param()

    call read_3body_param()

    call returnclose_iounit(ff_unit)
    call returnclose_iounit(input_epe)

    return

100 call wrong_line(epeinp)

  contains
!....................................................................................
    subroutine get_type()

      integer(i4_kind) :: ig,ig1,status,iat,ntp
      logical :: previous_type

      ntp=size(type_nm)
      n_type_atm=0
      
      allocate(at_type(n_atm),stat=status)
      ASSERT(status==0)

      do ig=1,n_atm
         if(n_type_atm == 0) then
            at_type(ig)=1
            n_type_atm=1
         else
            previous_type=.false.
            do ig1=1,ig-1
               if(atnm(ig)==atnm(ig1)) then
                  if(ig > n_cluster .and. ig1 <= n_cluster) call error_handler( &
                       "Read_interface_param:Coincided atom names in GXFILE and EPE.PCS files")
                  previous_type=.true.
                  at_type(ig)=at_type(ig1)
                  exit
               end if
            end do
            if(.not.previous_type) then
               n_type_atm=n_type_atm+1
               at_type(ig)=n_type_atm
            end if
         end if
         iat=at_type(ig)
         ASSERT(iat <= ntp)
         type_nm(iat)=atnm(ig)
      end do

    end subroutine get_type
!....................................................................................
    subroutine read_2body_param()

      character(len=80) :: buffer
      character(len=6) :: atm_name(2)
      real(r8_kind) :: b0,r0,c0,cutoff,r_at_name(2)
      real(r8_kind), parameter :: small=1.0e-3_r8_kind
      integer(i4_kind) :: ir,ir1,jr,jr1,n_2b,status

      allocate(interface_ff(n_type_atm,n_type_atm),stat=status)
      ASSERT(status==0)

      do ir=1,n_type_atm
         do ir1=1,n_type_atm
            interface_ff(ir,ir1)%b=0.0_r8_kind
            interface_ff(ir,ir1)%r=1.0_r8_kind
            interface_ff(ir,ir1)%c=0.0_r8_kind
            interface_ff(ir,ir1)%d=0.0_r8_kind
            interface_ff(ir,ir1)%cut=0.0_r8_kind     
         end do
      end do

      read_block=find_block(ff_unit,"&TWO-BODY",fflib)
      if(.not.read_block) then
         call error_handler("Read_2body_param: No atomic two body parameters found")
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

            do jr=1,n_type_atm
               if(abs(r_at_name(1)- type_nm(jr)) <= small) exit
            end do
            do jr1=1,n_type_atm
               if(abs(r_at_name(2)-type_nm(jr1)) <= small) exit
            end do
            if(jr > n_type_atm .or. jr1 > n_type_atm) cycle
            n_2b=n_2b+1

            interface_ff(jr,jr1)%b=b0
            interface_ff(jr1,jr)%b=b0
            interface_ff(jr,jr1)%r=r0
            interface_ff(jr1,jr)%r=r0
            interface_ff(jr,jr1)%c=c0
            interface_ff(jr1,jr)%c=c0
            interface_ff(jr,jr1)%d=0.0_r8_kind
            interface_ff(jr1,jr)%d=0.0_r8_kind
            interface_ff(jr,jr1)%cut=cutoff    
            interface_ff(jr1,jr)%cut=cutoff    

         end do
      end if

      if(n_2b==0) then
         call error_handler( &
              "Read_2body_param: No 2body interaction between QM cluster and EPE")
      else
         write(*,*) '----------------------------------------------------------------------'
         write(*,*) '      The Potential Parameters acting between QM cluster and EPE'
         write(*,*) '----------------------------------------------------------------------'
         write(*,*) '                        Two-body parameters'
         write(*,*) '----------------------------------------------------------------------'
         write(*,*) '                   b           ro          c           d        cutoff'
         do ir=1,n_type_atm
            do jr=ir,n_type_atm
               if(interface_ff(ir,jr)%cut == 0.0_r8_kind) cycle 
               write(*,'(2f6.2,5f12.5)') type_nm(ir),type_nm(jr), &
                    interface_ff(ir,jr)%b,interface_ff(ir,jr)%r,interface_ff(ir,jr)%c, &
                    interface_ff(ir,jr)%d,interface_ff(ir,jr)%cut
            end do
         enddo
         write(*,*) '----------------------------------------------------------------------'
      end if

    return

100 call wrong_line(fflib)

    end subroutine read_2body_param
!....................................................................................

    subroutine read_3body_param()

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

      n_types_central_atoms_3body=0

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
                  do k=1,n_type_atm
                     if(abs(r_c_nm-type_nm(k)) <= small) then 
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
                  allocate(ki(n_type_atm,n_type_atm,n_type_atm), &
                       theta_0(n_type_atm,n_type_atm,n_type_atm),stat=status)
                  ASSERT(status.eq.0)
                  allocate(types(n_types_central_atoms_3body,5), r3b(n_types_central_atoms_3body), &
                       n_3b(n_types_central_atoms_3body),indexx(n_types_central_atoms_3body,20,3), &
                       stat=status)
                  ASSERT(status.eq.0)
                  types=0
                  ki=0.0_r8_kind
                  n_3b = 0
                  types = 0
               end if

               read(buffer,*,err=100) atm_nm,k_i,theta_i,r_3b
               read(atm_nm(1),*,iostat=status) at_name(1)
               if(status /= 0) cycle cycle1
               read(atm_nm(2),*,iostat=status) at_name(2)
               if(status /= 0) cycle cycle1
               read(atm_nm(3),*,iostat=status) at_name(3)
               if(status /= 0) cycle cycle1

               i1=0; i2=0; i3=0
               do k=1,n_type_atm
                  if(abs(at_name(1) - type_nm(k)) <= small) i1=k
                  if(abs(at_name(2) - type_nm(k)) <= small) i2=k
                  if(abs(at_name(3) - type_nm(k)) <= small) i3=k
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
         end do cycle1
      end if

      if (n_types_central_atoms_3body /= 0) then
         write(*,*) '              Three-body parameters'
         write(*,*) '------------------------------------------------------'
         write(*,*) '                       ki        theta         r3b'
         do i=1,n_types_central_atoms_3body
            do j=1,n_3b(i)
               write(*,'(3f6.2,3f12.5)') type_nm(indexx(i,j,1)), &
                    type_nm(indexx(i,j,2)),type_nm(indexx(i,j,3)), &
                    ki(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)), &
                    theta_0(indexx(i,j,1),indexx(i,j,2),indexx(i,j,3)), r3b(i)
            enddo
         enddo
         deallocate(n_3b,indexx,stat=status)
         ASSERT(status == 0)
      endif
      write(*,*) '------------------------------------------------------'

      return

100   call wrong_line(fflib)

    end subroutine read_3body_param
  end subroutine read_interface_param
!====================================================================================

!====================================================================================
  subroutine free_ff_arrays

    integer(i4_kind) :: status

    if(allocated(at_type)) then
       deallocate(at_type,stat=status)
       ASSERT(status==0)
    endif

    if(allocated(interface_ff)) then
       deallocate(interface_ff,stat=status)
       ASSERT(status==0)
    endif

  end subroutine free_ff_arrays
end module qm_epe_interface_module

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
program def_tip_waters

  implicit none

  character*256 :: buffer
  integer, parameter :: N_atoms=1000
  character*4 :: atname(N_atoms)
  real*8 :: r(3,N_atoms)
  character*2 :: w_name(4),w_name_s(5)
  integer*4 :: z_mat(4,3),z_mat_s(5,3)
  real*8 :: z_mat_dat(4,3),z_mat_dat_s(5,3)
  real*8, parameter :: tip_bond=0.9572d0
  real*8, parameter :: tip_angle=104.52d0
  real*8, parameter :: tip4p_bond=0.15d0
  logical :: write_gx,write_xyz
  character*80 :: xyz_name,gx_name
  character*6, parameter :: default_name='output'
  integer*4 :: n_at,line_num
  integer*4 :: z_out_mat(N_atoms,3),n_out_mat(N_atoms,3)

  real*8, parameter :: pi=3.1415926535897932368d0
  real*8, parameter :: pi_degree=180.0d0
  real*8, parameter :: deg2rad=pi/pi_degree
  real*8, parameter :: rad2deg=pi_degree/pi
  real*8, parameter :: b2a=0.529177249d0
  real*8, parameter :: a2b=1.889725988579d0

  integer*4, parameter :: n_names=99
  character*2 :: name(n_names)= &
       (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne", &
         "Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca", &
         "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
         "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr", &
         "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn", &
         "Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd", &
         "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tu","Yb", &
         "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg", &
         "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", &
         "Pa","U ","Np","Pu","Am","Cm","Bk","Cf","xx"/)
  character*2 :: name1(n_names)= &
       (/"H ","HE","LI","BE","B ","C ","N ","O ","F ","NE", &
         "NA","MG","AL","SI","P ","S ","CL","AR","K ","CA", &
         "SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN", &
         "GA","GE","AS","SE","BR","KR","RB","SR","Y ","ZR", &
         "NB","MO","TC","RU","RH","PD","AG","CD","IN","SN", &
         "SB","TE","I ","XE","CS","BA","LA","CE","PR","ND", &
         "PM","SM","EU","GD","TB","DY","HO","ER","TU","YB", &
         "LU","HF","TA","W ","RE","OS","IR","PT","AU","HG", &
         "TL","PB","BI","PO","AT","RN","FR","RA","AC","TH", &
         "PA","U ","NP","PU","AM","CM","BK","CF","XX"/)
  integer*4 :: icl,its,it,iou

  write_gx=.false.
  write_xyz=.true.

  icl=0; its=4000; it=4000; iou=4000 
  z_out_mat=0
  n_out_mat=0

  n_at=0
  line_num=0
  do 
     buffer=''
     if(eof(5)) goto 100
     line_num=line_num+1
     read(5,'(a256)',err=200) buffer
     if(get_word(buffer,"&CLUSTER")) then
        icl=line_num
        if(icl > its .or. icl > it .or. icl > iou) &
             call error(icl,'&CLUSTER has to be the first section in input file')
        call read_cluster()
     else if(get_word(buffer,"&TIP_SYM")) then
        its=line_num
        call read_water_sym()
        call add_water_to_cluster('sym')
     else if(get_word(buffer,"&TIP")) then
        it=line_num
        call read_water()
        call add_water_to_cluster('nosym')
     else if(get_word(buffer,"&OUTPUT")) then
        iou=line_num
        call read_output()
     else if(icl==0 .and. its==4000 .and. it==4000 .and. iou==4000) then
        call error(line_num,'Input file does not contain data')
     end if
  end do

200 call error(line_num,'Wrong input line')

100 continue

  if(line_num==0) call error(line_num,'Input file does not contain data')

  if(write_xyz) call xyz_output()
  if(write_gx) call gx_output()

contains

  subroutine error(inp_line,message)

    integer*4 :: inp_line
    character(len=*) ::message
    character*4 :: buf

    write(buf,'(i4)') inp_line

    write(*,*) message//'. '//buf//'-input line.'

    stop

  end subroutine error

  subroutine upcase(string)

    character(len=*) :: string
    integer*4 :: i,ln,ich

    ln = len(trim(string))
    do i=1,ln
       ich=iachar(string(i:i))
       if(ich>=97 .and. ich<=122) then
          ich=ich-32
          string(i:i)=achar(ich)
       end if
    end do

  end subroutine upcase

  function get_word(buf,word)
    logical :: get_word

    character(len=*) :: buf,word

    call upcase(buf)

    if(index(buf,trim(word)) /= 0) then
       get_word=.true.
    else
       get_word=.false.
    end if
    
  end function get_word

  subroutine read_cluster()

    character*256 :: buf
    real*8 :: convert
    integer*4 :: ir

    convert=1.0d0

    do
       line_num=line_num+1
       read(5,'(a256)',err=200) buf
       if(get_word(buf,"&END")) then 
          exit
       else if(get_word(buf,"A.U.")) then
          convert=b2a
       else if(eof(5)) then
          call error(line_num,'Wrong section &CLUSTER')
       else
          n_at=n_at+1
          read(buf,*,err=200) atname(n_at),r(:,n_at)
       end if
    end do

    do ir=1,n_at
       r(:,ir)=r(:,ir)*convert
    end do

    return

200 call error(line_num,'Wrong section &CLUSTER')

  end subroutine read_cluster

  subroutine read_water_sym()

    character*256 :: buf
    real*8, parameter :: small=1.0d-10
    integer*4 :: ii

    ii=0
    do
       line_num=line_num+1
       read(5,'(a256)',err=200) buf
       if(get_word(buf,"&END")) then 
          exit
       else if(eof(5)) then
          call error(line_num,'&TIP_SYMR section is not terminated (&END)')
       else
          ii=ii+1
          call upcase(buf)
          if(ii==1) then 
             read(buf,*,err=200) w_name_s(1),z_mat_s(1,1:3),z_mat_dat_s(1,1:3)
          else if(ii==2) then
             read(buf,*,err=200) w_name_s(2),z_mat_s(2,2:3),z_mat_dat_s(2,2:3)
             z_mat_s(2,1)=n_at+1; z_mat_dat_s(2,1)=tip4p_bond
             w_name_s(3)='XX'
             z_mat_s(3,1)=n_at+2; z_mat_s(3,2)=n_at+1; z_mat_s(3,3)=z_mat_s(1,1)
             z_mat_dat_s(3,1)=0.5d0; z_mat_dat_s(3,2)=90.0d0; z_mat_dat_s(3,3)=0.0d0
             w_name_s(4)='H'
             z_mat_s(4,1)=n_at+1; z_mat_s(4,2)=n_at+2; z_mat_s(4,3)=n_at+3
             z_mat_dat_s(4,1)=tip_bond; z_mat_dat_s(4,2)=tip_angle*0.5d0; z_mat_dat_s(4,3)=90.0d0
             w_name_s(5)='H'
             z_mat_s(5,1)=n_at+1; z_mat_s(5,2)=n_at+2; z_mat_s(5,3)=n_at+3
             z_mat_dat_s(5,1)=tip_bond; z_mat_dat_s(5,2)=tip_angle*0.5d0; z_mat_dat_s(5,3)=-90.0d0
          end if
       end if
    end do

    return

200 call error(line_num,'Wrong section &TIP_SYM')

  end subroutine read_water_sym

  subroutine read_water()

    character*256 :: buf
    real*8, parameter :: small=1.0d-10
    integer*4 :: ii

    ii=0
    do
       line_num=line_num+1
       read(5,'(a256)',err=200) buf
       if(get_word(buf,"&END")) then 
          exit
       else if(eof(5)) then
          call error(line_num,'&TIP section is not terminated (&END)')
       else
          ii=ii+1
          call upcase(buf)
          if(ii==1) then
             read(buf,*,err=200) w_name(1),z_mat(1,1:3),z_mat_dat(1,1:3)
          else if(ii==2) then
             read(buf,*,err=200) w_name(2),z_mat(2,2:3),z_mat_dat(2,2:3)
             z_mat(2,1)=n_at+1; z_mat_dat(2,1)=tip_bond
          else if(ii==3) then
             read(buf,*,err=200) w_name(4),z_mat(4,3),z_mat_dat(4,3)
             z_mat(4,1)=n_at+2; z_mat(4,2)=n_at+1
             z_mat_dat(4,1)=tip_bond; z_mat_dat(4,2)=tip_angle
             w_name(3)='OQ'
             z_mat(3,1)=n_at+2; z_mat(3,2)=n_at+1; z_mat(3,3)=z_mat(4,3)
             z_mat_dat(3,1)=tip4p_bond; z_mat_dat(3,2)=tip_angle*0.5d0;
             z_mat_dat(3,3)=z_mat_dat(4,3)
          end if
       end if
    end do

    return

200 call error(line_num,'Wrong section &TIP')

  end subroutine read_water

  subroutine add_water_to_cluster(w_type)

    character(len=*) :: w_type
    real*8 :: xyz(4,3)
    integer*4 :: a1,a2,a3
    real*8 :: rw(3,3)
    integer*4 :: ia
    
    if(trim(w_type) == 'sym') then
       do ia=1,5
          a3=z_mat_s(ia,1); a2=z_mat_s(ia,2); a1=z_mat_s(ia,3)
          rw(:,1)=r(:,a1); rw(:,2)=r(:,a2); rw(:,3)=r(:,a3)
          call add_1_4_atom(xyz(ia,:),rw,z_mat_dat_s(ia,:))
          n_at=n_at+1
          r(:,n_at)=xyz(ia,:)
          atname(n_at)=w_name_s(ia)
          z_out_mat(n_at,:)=z_mat_s(ia,:)
          if(ia==1) then
             n_out_mat(n_at,:)=1
          else if(ia==2) then
             n_out_mat(n_at,2:3)=1
          else if(ia==3) then
             n_out_mat(n_at,3)=1
          end if
       end do
    else if((w_type) == 'nosym') then
       do ia=1,4
          a3=z_mat(ia,1); a2=z_mat(ia,2); a1=z_mat(ia,3)
          rw(:,1)=r(:,a1); rw(:,2)=r(:,a2); rw(:,3)=r(:,a3)
          call add_1_4_atom(xyz(ia,:),rw,z_mat_dat(ia,:))
          n_at=n_at+1
          r(:,n_at)=xyz(ia,:)
          atname(n_at)=w_name(ia)
          if(ia/=4) then
             z_out_mat(n_at,:)=z_mat(ia,:)
          else
             z_out_mat(n_at,1)=n_at-1; z_out_mat(n_at,2:3)=z_mat(ia-1,1:2) 
          end if
          if(ia==1) then
             n_out_mat(n_at,:)=1
          else if(ia==2) then
             n_out_mat(n_at,2:3)=1
          else if(ia==3) then
             n_out_mat(n_at,3)=1
          end if
       end do
    end if

  end subroutine add_water_to_cluster

  subroutine add_1_4_atom(ra,rb,internals)

    real*8 :: ra(3),rb(3,3)
    real*8 :: internals(3)
    real*8 :: bond,angle,dihedral
    real*8 :: cos_a,sin_a,cos_d,sin_d
    real*8 :: r21(3),d21,r32(3),d32,r1(3),d1,r2(3),d2,r3(3),r4(3),d4
    real*8, parameter :: small=1.0d-10

    bond=internals(1)
    angle=internals(2)*deg2rad
    dihedral=internals(3)*deg2rad

    cos_a=cos(angle); sin_a=sin(angle)
    cos_d=cos(dihedral); sin_d=sin(dihedral)

    r21=rb(:,2)-rb(:,1); r32=rb(:,3)-rb(:,2)
    d21=sqrt(dot_product(r21,r21))
    r21=r21/d21
    d32=sqrt(dot_product(r32,r32))
    r32=r32/d32

    r1=vector_product(r21,r32)
    d1=sqrt(dot_product(r1,r1))
    r1=r1/d1

    r2=vector_product(r1,r32)
    d2=sqrt(dot_product(r2,r2))
    r2=r2/d2

    r3=r2*cos_d+r1*sin_d

    r4=-r32*cos_a+r3*sin_a
    d4=sqrt(dot_product(r4,r4))
    r4=r4*bond/d4

    ra=rb(:,3)+r4

  end subroutine add_1_4_atom

  function vector_product(v1,v2)

    real*8 :: vector_product(3)
    real*8 :: v1(3),v2(3)

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product

  subroutine read_output()

    character*256 :: buf,buf1
    character*3 :: abc

    do
       line_num=line_num+1
       read(5,'(a256)',err=300) buf
       buf1=buf
       if(get_word(buf,"&END")) then 
          exit
       else if(get_word(buf,"XYZ")) then
          write_xyz=.true.
          read(buf1,*,err=100) abc,xyz_name
          xyz_name=trim(xyz_name)//'.xyz'
          cycle
100       xyz_name=trim(default_name)//'.xyz'
       else if(get_word(buf,"GX")) then
          write_gx=.true.
          read(buf1,*,err=200) abc,gx_name
          gx_name='gx.'//trim(gx_name)
          cycle
200       gx_name='gx.'//trim(default_name)
       else if(eof(5)) then
          call error(line_num,'&OUTPUT section is not terminated (&END)')
       end if
    end do

    return

300 call error(line_num,'Wrong section &OUTPUT')

  end subroutine read_output

  subroutine xyz_output()

    integer*4 :: io
    character*4 :: an

    open(10,file=trim(xyz_name))
    write(10,*) n_at
    write(10,*) '!!!!!!!!!!!!!!'
    do io=1,n_at
       an=atname(io)
       if(trim(an)=='OQ') an='O'
       write(10,'(a4,3f24.12)') an,r(:,io)
    end do    
    close(10)

  end subroutine xyz_output

  subroutine gx_output()

    integer*4 :: io,jo
    character*2 :: an
    real*8 :: at_num
    integer*4 :: num,numa,num1

     open(10,file=trim(gx_name))
     num=0; num1=0
     do io=1,n_at
        an=atname(io)
        if(an=="OQ") then
           at_num=8.01
        else
           do jo=1,n_names
              if(an==name(jo) .or. trim(an)==name1(jo)) then
                 at_num=float(jo)
                 exit
              end if
           end do
        end if
        num1=num1+1
        if(at_num /= 99.0d0) then
           num=num+1
           numa=num
        else
           numa=0
        end if
        write(10,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4)') &
             at_num,r(:,io)*a2b,numa,num1,z_out_mat(io,:),n_out_mat(io,:)
     end do
     write(10,'(f5.2)') -1.0
     close(10)

  end subroutine gx_output

end program def_tip_waters

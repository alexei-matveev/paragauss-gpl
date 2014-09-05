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
program def_tip_from_xyz

  implicit none

  character*80 :: xyz_file1, xyz_file2
  integer*4 :: n_cluster, n_total
  integer, parameter :: N_atoms=1000
  character*4 :: atname(N_atoms),atname1(N_atoms)
  real*8 :: r(3,N_atoms),r1(3,N_atoms) 
  integer*4 :: z_out_mat(N_atoms,3),n_out_mat(N_atoms,3)

  real*8, parameter :: tip_bond=0.9572d0
  real*8, parameter :: tip_angle=104.52d0
  real*8, parameter :: tip4p_bond=0.15d0
  real*8, parameter :: water_bond=1.1d0

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

  character*40 :: buf
  logical :: yes
  integer*4 :: i,j,k

  z_out_mat=0; n_out_mat=0
  r=0.0d0; r1=0.0d0

  call getarg(1,xyz_file1)
  call getarg(2,xyz_file2)

  if(trim(xyz_file1)=='' .or. trim(xyz_file2)=='') then
     write(*,*) 'Please, give names of XYZ files'
     stop
  end if

  inquire(file=trim(xyz_file1),exist=yes)
  if(.not.yes) then
     write(*,*) 'XYZ file '//trim(xyz_file1)//' does not exist'
     stop
  end if
  inquire(file=trim(xyz_file2),exist=yes)
  if(.not.yes) then
     write(*,*) 'XYZ file '//trim(xyz_file2)//' does not exist'
     stop
  end if

  !reading in atomic coordinates
  open(10,file=trim(xyz_file1)); open(11,file=trim(xyz_file2))

  read(10,*) n_cluster
  read(11,*) n_total

  if(n_cluster >= n_total) then
     write(*,*) 'The first XYZ file is larger than the second one'
     stop
  end if

  read(11,'(a40)') buf
100 do i=1,n_total
     read(11,*) atname(i),r(:,i)
     call upcase(atname(i))
  end do
  close(11)

  read(10,'(a40)') buf
200 do i=1,n_cluster
     read(10,*) atname(i),r(:,i)
     call upcase(atname(i))
  end do
  close(10)

  call check_waters()

  !adding tip4p atoms
  i=0; j=0; k=0
  do
     i=i+1; j=j+1
     if(i > n_total) exit
     if(i<=n_cluster) then
        atname1(j)=atname(i)
        r1(:,j)=r(:,i)
     else
        k=k+1
        if(k==1 .and. trim(atname(i))=='O') then
           atname1(j)=atname(i)
           r1(:,j)=r(:,i)
           j=j+1
           atname1(j)='OQ'
           r1(:,j)=0.0d0
           j=j+1
           atname1(j)='XX'
           r1(:,j)=0.0d0
        else if(k==2 .and. trim(atname(i))=='O') then
           atname1(j)=atname(i)
           r1(:,j)=r(:,i)
           j=j+1
           atname1(j)='OQ'
           r1(:,j)=0.0d0
        else
           atname1(j)=atname(i)
           r1(:,j)=r(:,i)
        end if
        if(k==3) k=0
     end if
  end do
  n_total=j-1

  call make_tip4p()

  call gx_output()

contains

  subroutine upcase(string)

    character(len=*) :: string
    integer*4 :: iu,ln,ich

    ln = len(trim(string))
    do iu=1,ln
       ich=iachar(string(iu:iu))
       if(ich>=97 .and. ich<=122) then
          ich=ich-32
          string(iu:iu)=achar(ich)
       end if
    end do

  end subroutine upcase

  subroutine check_waters()

    real*8 :: r12(3),d12,r32(3),d32
    integer*4 :: i1,i2,i3
    character*4 :: buf
    integer*4 :: ic,jc,kc

    jc=0
    do ic=n_cluster+1,n_total,3
       jc=jc+1
       write(buf,'(i4)') jc
       do kc=1,3
          if(trim(atname(ic+kc-1)) /= 'O' .and. &
               trim(atname(ic+kc-1)) /= 'H') goto 100
       end do
       if(trim(atname(ic)) == 'H' .and. trim(atname(ic+1)) /= 'O') goto 100
       if(trim(atname(ic)) == 'O' .and. trim(atname(ic+1)) /= 'H') goto 100
       if(trim(atname(ic+2)) == 'O') goto 100

       if(trim(atname(ic)) == 'O') then
          i1=ic+1; i2=ic; i3=ic+2
       else if(trim(atname(ic+1)) == 'O') then
          i1=ic; i2=ic+1; i3=ic+2
       end if
       r12=r(:,i1)-r(:,i2); d12=sqrt(dot_product(r12,r12))
       r32=r(:,i3)-r(:,i2); d32=sqrt(dot_product(r32,r32))
       if(d12 >= water_bond .or.d32 >= water_bond) then 
          goto 200
       else
          if(d12 /= tip_bond) then
             r12=r12*tip_bond/d12
             r(:,i1)=r(:,i2)+r12
          end if
          if(d32 /= tip_bond) then
             r32=r32*tip_bond/d32
             r(:,i3)=r(:,i2)+r32
          end if
       end if
    end do

    return

100 write(*,*) 'The '//buf//'th water are not water'
    stop
200 write(*,*) 'The '//buf//'th water: Atoms are too far from each other'
    write(*,*) 'Is this water?'
    stop

  end subroutine check_waters

  subroutine make_tip4p()

    real*8 :: rq(3),rq2(3),dq2,rm(3,3)
    real*8 :: zmat(3)
    integer*4 :: im,nw

    im=n_cluster; nw=1
    do
       im=im+nw
       if(im > n_total) exit
       if(trim(atname1(im)) == 'O') then
          nw=5
          rq=(r1(:,im+3)+r1(:,im+4))/2.0d0
          rq2=rq-r1(:,im); dq2=sqrt(dot_product(rq2,rq2))
          rq2=rq2*tip4p_bond/dq2; r1(:,im+1)=r1(:,im)+rq2

          rm(:,1)=r1(:,im+3); rm(:,2)=r1(:,im); rm(:,3)=r1(:,im+1)
          zmat(1)=0.5d0;  zmat(2)=90.0d0; zmat(3)=90.0d0
          call add_1_4_atom(r1(:,im+2),rm,zmat)

          rm(:,1)=r1(:,im+2); rm(:,2)=r1(:,im+1); rm(:,3)=r1(:,im)
          zmat(1)=tip_bond;  zmat(2)=tip_angle*0.5d0; zmat(3)=90.0d0
          call add_1_4_atom(r1(:,im+3),rm,zmat)

          rm(:,1)=r1(:,im+2); rm(:,2)=r1(:,im+1); rm(:,3)=r1(:,im)
          zmat(1)=tip_bond;  zmat(2)=tip_angle*0.5d0; zmat(3)=-90.0d0
          call add_1_4_atom(r1(:,im+4),rm,zmat)

          n_out_mat(im,:)=1
          z_out_mat(im+1,1)=im; n_out_mat(im+1,2:3)=1
          z_out_mat(im+2,1)=im+1; z_out_mat(im+2,2)=im; n_out_mat(im+2,3)=1
          z_out_mat(im+3,1)=im; z_out_mat(im+3,2)=im+1; z_out_mat(im+3,3)=im+2
          z_out_mat(im+4,1)=im; z_out_mat(im+4,2)=im+1; z_out_mat(im+4,3)=im+2
       else if(trim(atname1(im)) == 'H') then
          nw=4
          rq=(r1(:,im)+r1(:,im+3))/2.0d0
          rq2=rq-r1(:,im+1); dq2=sqrt(dot_product(rq2,rq2))
          rq2=rq2*tip4p_bond/dq2; r1(:,im+2)=r1(:,im+1)+rq2

          rm(:,1)=r1(:,im+3); rm(:,2)=r1(:,im+2); rm(:,3)=r1(:,im+1)
          zmat(1)=tip_bond;  zmat(2)=tip_angle*0.5d0; zmat(3)=180.0d0
          call add_1_4_atom(r1(:,im),rm,zmat)

          rm(:,1)=r1(:,im); rm(:,2)=r1(:,im+2); rm(:,3)=r1(:,im+1)
          zmat(1)=tip_bond;  zmat(2)=tip_angle*0.5d0; zmat(3)=-180.0d0
          call add_1_4_atom(r1(:,im+3),rm,zmat)

          n_out_mat(im,:)=1
          z_out_mat(im+1,1)=im; n_out_mat(im+1,2:3)=1
          z_out_mat(im+2,1)=im+1; z_out_mat(im+2,2)=im; n_out_mat(im+2,3)=1
          z_out_mat(im+3,1)=im+2; z_out_mat(im+3,2)=im+1; z_out_mat(im+3,3)=im
       end if
    end do

  end subroutine make_tip4p

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

  subroutine gx_output()

    integer*4 :: io,jo
    character*2 :: an
    real*8 :: at_num
    integer*4 :: num,numa,num1

     num=0; num1=0
     do io=1,n_total
        an=atname1(io)
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
        write(6,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4)') &
             at_num,r1(:,io)*a2b,numa,num1,z_out_mat(io,:),n_out_mat(io,:)
     end do
     write(6,'(f5.2)') -1.0

  end subroutine gx_output

end program def_tip_from_xyz

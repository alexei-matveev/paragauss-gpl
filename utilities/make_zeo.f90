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
program make_zeo

  implicit none

  character*80 :: line
  character*37, parameter :: start_line = &
       "Final fractional coordinates of atoms"
  character*33, parameter :: start_line2 = &
       "Final asymmetric unit coordinates"
  character*10, parameter :: final_line = &
       "----------"
  character*37, parameter :: start_line1 = &
       "Final cell parameters and derivatives"
  character*10 :: buffer
  integer :: i,j
  integer :: counter_N,counter_F,counter_Al,counter_P,counter_B

  integer :: N_atoms,N_at,N_core
  real(8) :: r(3),dr,rf(3),ra(3),rp(3),rh(3)

  type atom
     integer :: num
     character*3 :: name
     character*1 :: c_s
     real(8) :: r(3)
  end type atom
  type(atom) :: atoms_fra(1000)
  type(atom) :: atoms_car(1000)

  type cell
     real(8) :: a(6)
  end type cell
  type(cell) :: latt

  type vectors
     real(8) :: v1(3),v2(3),v3(3)
  end type vectors
  type(vectors) :: l_vect

  real(8), parameter :: pi=3.141592653589793d0
  real(8), parameter :: pi_degree=180.0d0
  real(8), parameter :: Si_O=2.00d0

  character*255 :: gulp_out
  integer :: ff
  !==========================================================

  write(6,*) "Name of Gulp output file:"
  read(5,*) gulp_out
  write(6,*) "What does Force fiels have to be applied for?"
  write(6,*) "Either Sierka-Sauer DFT (1) or our FF3 extended (2):"
  read(5,*) ff

  if(ff /= 1 .and. ff /= 2) then
     print*,'Wrong parameter to chose Force Field ',ff
     stop
  end if

  open(15,file=trim(gulp_out))

  l1: do
     read(15,'(a80)') line
     if(index(line,start_line) /= 0 .or. &
          index(line,start_line2) /= 0) then
        l2: do i=1,5
           read(15,'(a80)') line
        end do l2
        exit l1
     end if
  end do l1
 
  i=0; j=0
  l3: do
     read(15,'(a80)') line
     if(index(line,final_line) /= 0) exit l3
     i=i+1
     read(line,*) atoms_fra(i)%num,atoms_fra(i)%name,atoms_fra(i)%c_s, &
          atoms_fra(i)%r
     if(atoms_fra(i)%c_s == "c") j=j+1
!!$     print*,atoms_fra(i)%num,atoms_fra(i)%name,atoms_fra(i)%c_s, &
!!$          atoms_fra(i)%r 
  end do l3
  N_atoms=i; N_core=j

  l4: do
     read(15,'(a80)') line
     if(index(line,start_line1) /= 0) then
        l5: do i=1,2
           read(15,'(a80)') line
        end do l5
        exit l4
     end if
  end do l4

  i=0
  l6: do
     read(15,'(a80)') line
     if(index(line,final_line) /= 0) exit l6
     i=i+1
     read(line,*) buffer,latt%a(i)
!!$     print*, latt%a(i)
  end do l6

  close(15)

  call cell2vect(latt,l_vect)

!!$  print*,l_vect%v1
!!$  print*,l_vect%v2
!!$  print*,l_vect%v3

  call fra2car(N_atoms,l_vect,atoms_fra,atoms_car)

!!$  do i=1,N_atoms
!!$     print*,atoms_car(i)%num,atoms_car(i)%name,atoms_car(i)%c_s, &
!!$          atoms_car(i)%r
!!$     if(atoms_car(i)%c_s == "s") cycle
!!$     print*,atoms_car(i)%name,atoms_car(i)%r
!!$  end do

  if(ff == 1) then
     call ssd()
  else
     call extend_ff3()
  end if

  atoms_fra(:)%name=atoms_car(:)%name
  atoms_fra(:)%c_s=atoms_car(:)%c_s

  do i=N_atoms+1,N_at
     call car2fra(l_vect,atoms_car(i)%r,atoms_fra(i)%r)
  end do

  call write_vm_car()
  call write_vm_fra()
  call write_xyz_car()
  if(ff == 1) then
     call gulp_inp_ssd()
  else
     call gulp_inp_ext_ff3()
  end if

contains

  subroutine ssd()

    l1: do i=1,N_atoms
       if(atoms_car(i)%c_s == "s") cycle l1
       if(trim(atoms_car(i)%name) == "F".or. &
           trim(atoms_car(i)%name) == "O") then
          rf=atoms_car(i)%r
          counter_Al=0
          l2: do j=1,N_atoms
             if(j==i .or. atoms_car(j)%c_s == "s") cycle l2
             r=atoms_car(j)%r-atoms_car(i)%r
             call image(r)
             dr=sqrt(dot_product(r,r))
             if(dr <= Si_O) then
                if(trim(atoms_car(j)%name) == "Al") then
                   counter_Al=counter_Al+1
                end if
             end if
          end do l2
          if(counter_Al > 1) then
             print*,"Oxigen with number ",atoms_car(i)%num, &
                  " contacts with two Al"
             stop
          end if
       end if
    end do l1

    N_at=N_atoms
    l9: do i=1,N_atoms
       if(atoms_car(i)%c_s == "s") cycle l9
       if(trim(atoms_car(i)%name) == "F") then
          rf=atoms_car(i)%r
          counter_P=0; counter_Al=0
          l10: do j=1,N_atoms
             if(j==i .or. atoms_car(j)%c_s == "s") cycle l10
             r=atoms_car(j)%r-atoms_car(i)%r
             call image(r)
             dr=sqrt(dot_product(r,r))
             if(dr <= Si_O) then
                if(trim(atoms_car(j)%name) == "Al") then
!!$                   ra=atoms_car(j)%r
                   ra=rf+r
                   counter_Al=counter_Al+1
                elseif(trim(atoms_car(j)%name) == "Si") then
!!$                   rp=atoms_car(j)%r
                   rp=rf+r
                   counter_P=counter_P+1
                elseif(trim(atoms_car(j)%name) == "H") then
                   cycle l9 
                end if
             end if
          end do l10
          if(counter_Al /= 1 .and. counter_P /= 1) then
             print*,"Check the environment of O3(F) atom with number ", &
                  atoms_car(i)%num
             stop
          end if
          call get_H_position(ra,rf,rp,rh)
          N_at=N_at+1
          N_core=N_core+1
          atoms_car(N_at)%num=N_at
          atoms_car(N_at)%name="H"
          atoms_car(N_at)%c_s="c"
          atoms_car(N_at)%r=rh
       end if
    end do l9

  end subroutine ssd

  subroutine extend_ff3()

    !Search of Al, O3(F) and O2(N) atoms  
    counter_Al=0
    l7: do i=1,N_atoms
       if(atoms_car(i)%c_s == "s") cycle l7
       if(trim(atoms_car(i)%name) == "Al") then
          counter_Al=counter_Al+1
          counter_N=0; counter_F=0
          l8: do j=1,N_atoms
             if(j==i .or. atoms_car(j)%c_s == "s") cycle l8
             r=atoms_car(j)%r-atoms_car(i)%r
!!$           print*,i,j,sqrt(dot_product(r,r))
             call image(r)
             dr=sqrt(dot_product(r,r))
!!$           print*,dr
             if(dr <= Si_O) then
                counter_N=counter_N+1
                if(trim(atoms_car(j)%name) == "F") then
                   counter_F=counter_F+1
                   cycle l8
                elseif(trim(atoms_car(j)%name) == "N") then
                   print*,"Al atom with number ",atoms_car(i)%num, &
                        " located too close to previously defined Al atom"
                   stop
                elseif(trim(atoms_car(j)%name) == "O") then
                   atoms_car(j)%name="N"
                end if
             end if
!!$           print*,counter_N
          end do l8
          if(counter_F == 0) then
             print*,"Al atom with number ",atoms_car(i)%num, &
                  " has no a bond with O3(F) atom"
             stop
          elseif(counter_F > 1) then
             print*,"Al atom with number ",atoms_car(i)%num, &
                  " has ",counter_F," bonds with O3(F) atoms"
             stop
          end if
          if(counter_N /= 4) then
             print*,"Al atom with number ",atoms_car(i)%num, &
                  " has wrong oxygen environment"
             stop
          end if
       end if
    end do l7
    if(counter_Al == 0) then
       print*,"The system has no Al atoms. Nothing to do"
       stop
    end if

    !Search of Si2(P) and definition of H positions
    N_at=N_atoms
    l9: do i=1,N_atoms
       if(atoms_car(i)%c_s == "s") cycle l9
       if(trim(atoms_car(i)%name) == "F") then
          rf=atoms_car(i)%r
          counter_P=0; counter_Al=0
          l10: do j=1,N_atoms
             if(j==i .or. atoms_car(j)%c_s == "s") cycle l10
             r=atoms_car(j)%r-atoms_car(i)%r
             call image(r)
             dr=sqrt(dot_product(r,r))
             if(dr <= Si_O) then
                if(trim(atoms_car(j)%name) == "Al") then
!!$                   ra=atoms_car(j)%r
                   ra=rf+r
                   counter_Al=counter_Al+1
                else if(trim(atoms_car(j)%name) == "Si" .or. &
                     trim(atoms_car(j)%name) == "P") then
!!$                   rp=atoms_car(j)%r
                   rp=rf+r
                   counter_P=counter_P+1
                   atoms_car(j)%name="P"
                else if(trim(atoms_car(j)%name) == "H") then
                    cycle l9
                end if
             end if
          end do l10
          if(counter_Al /= 1 .and. counter_P /= 1) then
             print*,"Check the environment of O3(F) atom with number ", &
                  atoms_car(i)%num
             stop
          end if
          call get_H_position(ra,rf,rp,rh)
          N_at=N_at+1
          N_core=N_core+1
          atoms_car(N_at)%num=N_at
          atoms_car(N_at)%name="H"
          atoms_car(N_at)%c_s="c"
          atoms_car(N_at)%r=rh
       end if
    end do l9

    !Search O4(B) atoms
    l11: do i=1,N_atoms
       if(atoms_car(i)%c_s == "s") cycle l11
       if(trim(atoms_car(i)%name) == "P") then
          counter_B=0; counter_F=0;
          l12: do j=1,N_atoms
             if(j==i .or. atoms_car(j)%c_s == "s") cycle l12
             r=atoms_car(j)%r-atoms_car(i)%r
             call image(r)
             dr=sqrt(dot_product(r,r))
             if(dr <= Si_O) then
                if(trim(atoms_car(j)%name) == "F") then
                   counter_F=counter_F+1
                elseif(trim(atoms_car(j)%name) == "O") then
                   counter_B=counter_B+1
                   atoms_car(j)%name="B"
                elseif(trim(atoms_car(j)%name) == "B") then
                   counter_B=counter_B+1
                elseif(trim(atoms_car(j)%name) == "N") then
                   counter_B=counter_B+1
                end if
             end if
          end do l12
          if(counter_B+counter_F /= 4) then
             print*,"Check the environment of Si2(P) atom with number ", &
                  atoms_car(i)%num,counter_B+counter_F
             stop
          end if
       end if
    end do l11

  end subroutine extend_ff3

  subroutine cell2vect(ucell,vect)

    type(cell) :: ucell
    type(vectors) :: vect

    real(8) :: alpha_r, beta_r, gamma_r

    alpha_r=pi*ucell%a(4)/pi_degree
    beta_r=pi*ucell%a(5)/pi_degree
    gamma_r=pi*ucell%a(6)/pi_degree

    vect%v1(1)=ucell%a(1)
    vect%v1(2)=0.0d0
    vect%v1(3)=0.0d0
    vect%v2(1)=ucell%a(2)*cos(gamma_r)
    vect%v2(2)=ucell%a(2)*sin(gamma_r)
    vect%v2(3)=0.0d0
    vect%v3(1)=ucell%a(3)*cos(beta_r)
    vect%v3(2)=(ucell%a(2)*ucell%a(3)*cos(alpha_r)-vect%v2(1)*vect%v3(1))/vect%v2(2)
    vect%v3(3)=sqrt(ucell%a(3)*ucell%a(3)-vect%v3(1)*vect%v3(1)-vect%v3(2)*vect%v3(2))

  end subroutine cell2vect

  subroutine fra2car(natms,amat,frac,cart)

    integer :: natms
    type(vectors) :: amat
    type(atom) :: frac(natms)
    type(atom) :: cart(natms)

    real(8) :: bmat(3,3)
    integer :: ii

    cart=frac

    bmat(:,1)=amat%v1
    bmat(:,2)=amat%v2
    bmat(:,3)=amat%v3

    do ii=1,natms
       cart(ii)%r=matmul(bmat,frac(ii)%r)
    end do

  end subroutine fra2car

  subroutine car2fra(amat,cart,frac)

    type(vectors) :: amat
    real(8) :: cart(3)
    real(8)  :: frac(3)

    real(8) :: bmat(3,3),det,rr
    integer :: j1

    bmat(1,1)=amat%v2(2)*amat%v3(3)-amat%v2(3)*amat%v3(2)
    bmat(2,1)=amat%v1(3)*amat%v3(2)-amat%v1(2)*amat%v3(3)
    bmat(3,1)=amat%v1(2)*amat%v2(3)-amat%v1(3)*amat%v2(2)
    bmat(1,2)=amat%v2(3)*amat%v3(1)-amat%v2(1)*amat%v3(3)
    bmat(2,2)=amat%v1(1)*amat%v3(3)-amat%v1(3)*amat%v3(1)
    bmat(3,2)=amat%v1(3)*amat%v2(1)-amat%v1(1)*amat%v2(3)
    bmat(1,3)=amat%v2(1)*amat%v3(2)-amat%v2(2)*amat%v3(1)
    bmat(2,3)=amat%v1(2)*amat%v3(1)-amat%v1(1)*amat%v3(2)
    bmat(3,3)=amat%v1(1)*amat%v2(2)-amat%v1(2)*amat%v2(1)

    det=amat%v1(1)*bmat(1,1)+amat%v2(1)*bmat(2,1)+amat%v3(1)*bmat(3,1)
    rr=0.0d0
    if(abs(det) > 0.0d0) rr=1.0d0/det

    bmat=bmat*rr

    frac=matmul(bmat,cart)
!!$    do j1=1,3
!!$       if(abs(frac(j1)) /= 0.5d0) then
!!$          frac(j1)=frac(j1)-anint(frac(j1))
!!$       end if
!!$    end do

  end subroutine car2fra

  subroutine image(r1)

    real(8) :: r1(3)
    real(8) :: fr(3),bmat(3,3),det,rr
    integer :: ii

    bmat(1,1)=l_vect%v2(2)*l_vect%v3(3)-l_vect%v2(3)*l_vect%v3(2)
    bmat(2,1)=l_vect%v1(3)*l_vect%v3(2)-l_vect%v1(2)*l_vect%v3(3)
    bmat(3,1)=l_vect%v1(2)*l_vect%v2(3)-l_vect%v1(3)*l_vect%v2(2)
    bmat(1,2)=l_vect%v2(3)*l_vect%v3(1)-l_vect%v2(1)*l_vect%v3(3)
    bmat(2,2)=l_vect%v1(1)*l_vect%v3(3)-l_vect%v1(3)*l_vect%v3(1)
    bmat(3,2)=l_vect%v1(3)*l_vect%v2(1)-l_vect%v1(1)*l_vect%v2(3)
    bmat(1,3)=l_vect%v2(1)*l_vect%v3(2)-l_vect%v2(2)*l_vect%v3(1)
    bmat(2,3)=l_vect%v1(2)*l_vect%v3(1)-l_vect%v1(1)*l_vect%v3(2)
    bmat(3,3)=l_vect%v1(1)*l_vect%v2(2)-l_vect%v1(2)*l_vect%v2(1)

    det=l_vect%v1(1)*bmat(1,1)+l_vect%v2(1)*bmat(2,1)+l_vect%v3(1)*bmat(3,1)
    rr=0.0d0
    if(abs(det) > 0.0d0) rr=1.0d0/det

    bmat=bmat*rr

    fr=matmul(bmat,r1)
    do ii=1,3
       if(abs(fr(ii)) /= 0.5d0) then
          fr(ii)=fr(ii)-anint(fr(ii))
       end if
    end do

    bmat(:,1)=l_vect%v1
    bmat(:,2)=l_vect%v2
    bmat(:,3)=l_vect%v3

    r1=matmul(bmat,fr)

  end subroutine image

  subroutine get_H_position(r1,r2,r3,r4)

    real(8) :: r1(3),r2(3),r3(3),r4(3)
    real(8) :: r12(3),r32(3),dr12,dr32,drh
    real(8) :: cos_123,angle_123,angle_H,cos_H
    real(8) :: a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1
    real(8) :: m1,m2,m3,md
    real(8) :: maat(3,3)

    r12=r1-r2; r32=r3-r2
    dr12=sqrt(dot_product(r12,r12))
    dr32=sqrt(dot_product(r32,r32))
    cos_123=dot_product(r12,r32)/(dr12*dr32)
    angle_123=acos(cos_123)

    angle_H=(2.0d0*pi-angle_123)/2.0d0
    cos_H=cos(angle_H)

    drh=1.0d0

    maat(1,1)=r1(2); maat(1,2)=r1(3); maat(1,3)=1.0d0
    maat(2,1)=r2(2); maat(2,2)=r2(3); maat(2,3)=1.0d0
    maat(3,1)=r3(2); maat(3,2)=r3(3); maat(3,3)=1.0d0
    a1=det(maat)

    maat(1,1)=r1(3); maat(1,2)=r1(1); maat(1,3)=1.0d0
    maat(2,1)=r2(3); maat(2,2)=r2(1); maat(2,3)=1.0d0
    maat(3,1)=r3(3); maat(3,2)=r3(1); maat(3,3)=1.0d0
    b1=det(maat)

    maat(1,1)=r1(1); maat(1,2)=r1(2); maat(1,3)=1.0d0
    maat(2,1)=r2(1); maat(2,2)=r2(2); maat(2,3)=1.0d0
    maat(3,1)=r3(1); maat(3,2)=r3(2); maat(3,3)=1.0d0
    c1=det(maat)


    maat(1,1)=r1(1); maat(1,2)=r1(2); maat(1,3)=r1(3)
    maat(2,1)=r2(1); maat(2,2)=r2(2); maat(2,3)=r2(3)
    maat(3,1)=r3(1); maat(3,2)=r3(2); maat(3,3)=r3(3)
    d1=-det(maat)

    e1=r1(1)-r2(1); f1=r1(2)-r2(2); g1=r1(3)-r2(3)
    h1=-(r1(1)-r2(1))*r2(1)-(r1(2)-r2(2))*r2(2)-(r1(3)-r2(3))*r2(3)-dr12*drh*cos_H

    i1=r3(1)-r2(1); j1=r3(2)-r2(2); k1=r3(3)-r2(3)
    l1=-(r3(1)-r2(1))*r2(1)-(r3(2)-r2(2))*r2(2)-(r3(3)-r2(3))*r2(3)-dr32*drh*cos_H

    maat(1,1)=a1; maat(1,2)=b1; maat(1,3)=c1
    maat(2,1)=e1; maat(2,2)=f1; maat(2,3)=g1
    maat(3,1)=i1; maat(3,2)=j1; maat(3,3)=k1
    md=det(maat)
    
    maat(1,1)=-d1; maat(1,2)=b1; maat(1,3)=c1
    maat(2,1)=-h1; maat(2,2)=f1; maat(2,3)=g1
    maat(3,1)=-l1; maat(3,2)=j1; maat(3,3)=k1
    m1=det(maat)

    maat(1,1)=a1; maat(1,2)=-d1; maat(1,3)=c1
    maat(2,1)=e1; maat(2,2)=-h1; maat(2,3)=g1
    maat(3,1)=i1; maat(3,2)=-l1; maat(3,3)=k1
    m2=det(maat)

    maat(1,1)=a1; maat(1,2)=b1; maat(1,3)=-d1
    maat(2,1)=e1; maat(2,2)=f1; maat(2,3)=-h1
    maat(3,1)=i1; maat(3,2)=j1; maat(3,3)=-l1
    m3=det(maat)

    r4(1)=m1/md; r4(2)=m2/md; r4(3)=m3/md

  end subroutine get_H_position

  function det(mat)

    real(8) :: det
    real(8) :: mat(3,3)

    det= mat(1,1)*mat(2,2)*mat(3,3)+ &
         mat(1,2)*mat(2,3)*mat(3,1)+ &
         mat(1,3)*mat(2,1)*mat(3,2)- &
         mat(1,3)*mat(2,2)*mat(3,1)- &
         mat(1,2)*mat(2,1)*mat(3,3)- &
         mat(1,1)*mat(2,3)*mat(3,2)

  end function det

  subroutine write_vm_car()

    integer :: ij

    open(10,file="zeolite_c.dat")
 
    write(10,'(a7)')  " $title"
    write(10,'(a10)') " H-zeolite"
    write(10,'(a9)')  " $coord 1"
    do ij=1,N_at
       if(atoms_car(ij)%c_s == "s") cycle
       write(10,'(3f15.9,1x,a2)') atoms_car(ij)%r,trim(atoms_car(ij)%name)
    end do
    write(10,'(a11,6f10.5)') " $unitcell ",latt%a
    write(10,'(a5)')  " $end"
    
    close(10)

  end subroutine write_vm_car

  subroutine write_vm_fra()

    integer :: ij

    open(10,file="zeolite_f.dat")
 
    write(10,'(a7)')  " $title"
    write(10,'(a10)') " H-zeolite"
    write(10,'(a20)')  " $coord fractional 1"
    do ij=1,N_at
       if(atoms_fra(ij)%c_s == "s") cycle
       write(10,'(3f15.9,1x,a2)') atoms_fra(ij)%r,trim(atoms_fra(ij)%name)
    end do
    write(10,'(a11,6f10.5)') " $unitcell ",latt%a
    write(10,'(a5)')  " $end"
    
    close(10)

  end subroutine write_vm_fra

  subroutine write_xyz_car()

    integer :: ij

    open(10,file="zeolite.xyz")

    write(10,'(i5)') N_core
    write(10,'(a43)') "This file was generated by MAKE_ZEO program"
    do ij=1,N_at
       if(atoms_car(ij)%c_s == "s") cycle
       write(10,'(a2,3x,3f15.9)') trim(atoms_car(ij)%name),atoms_car(ij)%r
    end do

    close(10)

  end subroutine write_xyz_car

  subroutine gulp_inp_ext_ff3()

    integer :: ij,ji
    character*3 :: atname
    real(8) :: rcs(3), dcs
    real(8) :: c_s=0.25d0

    open(10,file="gulp_inp.ff3")

    write(10,'(a44)') "opti prop conp phon distance angles fix_molq"
    write(10,'(a5)') "title"
    write(10,'(a55)') "H-zeolite (correct the title if desirable) extended FF3"
    write(10,'(a3)') "end"
    write(10,'(a4)') "cell"
    write(10,'(6f11.6)') latt%a
    write(10,'(a10,i5)') "fractional", N_at

    do ij=1,N_at
       if(atoms_fra(ij)%c_s == "c") then
          select case(trim(atoms_fra(ij)%name))
          case ("Al")
             atname="Al1"
          case ("Si")
             atname="Si1"
          case ("F")
             atname="P3"
          case ("N")
             atname="P2"
          case ("P")
             atname="Cl2"
          case ("B")
             atname="O4"
          case ("O")
             atname="O1"
          case ("H")
             atname="H1"
          end select
          write(10,'(a3,1x,a4,1x,3F12.7)') atname, "core",atoms_fra(ij)%r
       elseif(atoms_fra(ij)%c_s == "s") then
          do ji=1,N_at
             if(atoms_fra(ji)%c_s == "s") cycle
             if(ji==ij) cycle
             rcs=atoms_car(ij)%r-atoms_car(ji)%r
             call image(rcs)
             dcs=sqrt(dot_product(rcs,rcs))
             if(dcs > c_s) cycle
             exit
          end do
          select case(trim(atoms_fra(ji)%name))
          case ("Al")
             atname="Al1"
          case ("Si")
             atname="Si1"
          case ("F")
             atname="P3"
          case ("N")
             atname="P2"
          case ("P")
             atname="Cl2"
          case ("B")
             atname="O4"
          case ("O")
             atname="O1"
          case ("H")
             atname="H1"
          end select
          write(10,'(a3,1x,a4,1x,3F12.7)') atname, "shel",atoms_fra(ij)%r
       end if
    end do

    write(10,'(a5)') "space"
    write(10,'(a1)') "1"
    write(10,'(a7)') "species"
    write(10,'(a13)') "Al1 core  1.1"
    write(10,'(a13)') "Si1 core  1.2"
    write(10,'(a13)') "Cl2 core  1.2"
    write(10,'(a15)') "O1  core  2.387"
    write(10,'(a15)') "O1  shel -2.987"
    write(10,'(a15)') "P2  core  2.691"
    write(10,'(a15)') "P2  shel -3.441"
    write(10,'(a15)') "P3  core  2.972"
    write(10,'(a15)') "P3  shel -3.372"
    write(10,'(a15)') "O4  core  2.387"
    write(10,'(a15)') "O4  shel -2.987"
    write(10,'(a14)') "H1  core  0.35"
    write(10,'(a1)') " "
    write(10,'(a8)') "rtol 1.0"
    write(10,'(a7)') "element"
    write(10,'(a16)') "covalent 13 0.01"
    write(10,'(a16)') "covalent  8 0.01"
    write(10,'(a16)') "covalent  1 0.01"
    write(10,'(a16)') "covalent 15 0.73"
    write(10,'(a16)') "covalent 17 1.70"
    write(10,'(a3)') "end"
    write(10,'(a1)') " "
    write(10,'(a4)') "buck"
    write(10,'(a57)') "O1  shel O1 shel 95169.354     0.199100 43.636    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "O1  shel O4 shel 95169.354     0.199100 43.636    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "O1  shel P2 shel 17471.593     0.268175 413.39    0.0 8.0"
    write(10,'(a10)') "buck intra"
    write(10,'(a57)') "O4  shel O4 shel 59194.280     0.238835 147.22    0.0 8.0"
    write(10,'(a10)') "buck inter"
    write(10,'(a57)') "O4  shel O4 shel 95169.354     0.199100 43.636    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "O1  shel P3 shel 38115.954     0.239040 125.58    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "O4  shel P3 shel 38115.954     0.239040 125.58    0.0 8.0"
    write(10,'(a10)') "buck inter"
    write(10,'(a57)') "P2  shel P2 shel 93947.963     0.223816 57.488    0.0 8.0"
    write(10,'(a10)') "buck intra"
    write(10,'(a57)') "P2  shel P2 shel 17471.593     0.268175 413.39    0.0 8.0"
    write(10,'(a10)') "buck inter"
    write(10,'(a57)') "O4  shel P2 shel 17471.593     0.268175 413.39    0.0 8.0"
    write(10,'(a10)') "buck intra"
    write(10,'(a57)') "P2  shel O4 shel 30471.574     0.259285 159.59    0.0 8.0"
    write(10,'(a10)') "buck intra"
    write(10,'(a57)') "P2  shel P3 shel 37214.853     0.238941 120.98    0.0 8.0"
    write(10,'(a10)') "buck inter"
    write(10,'(a57)') "P2  shel P3 shel 9418.5704     0.240302 6.4444    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "P3  shel H1 core 7518.5210     0.117834 2.0915    0.0 1.5"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Si1 core O1 shel 51431.799     0.174872 131.11    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Si1 core O4 shel 51431.799     0.174872 131.11    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Cl2 core O4 shel 55400.776     0.172102 144.27    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Si1 core P2 shel 35091.342     0.190375 174.57    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Cl2 core P2 shel 35291.445     0.189265 187.52    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Cl2 core P3 shel 40184.694     0.183955 166.99    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Al1 core P2 shel 29617.225     0.189495 122.39    0.0 8.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Al1 core P3 shel 40768.575     0.197235 181.04    0.0 8.0"
    write(10,'(a1)') " "
    write(10,'(a6)') "spring"
    write(10,'(a8)') "O1 78.05"
    write(10,'(a8)') "P2 84.76"
    write(10,'(a8)') "P3 88.45"
    write(10,'(a8)') "O4 78.05"
    write(10,'(a1)') " "
    write(10,'(a5)') "three"
    write(10,'(a53)') "P3  shel H1  core Cl2 core 9.7478  141.01 1.5 2.5 3.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "P3  shel H1  core Al1 core 9.6539  129.78 1.5 2.5 3.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Si1 core O1  shel O1  shel 0.89844 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Si1 core O1  shel O4  shel 0.89844 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Si1 core O1  shel P2  shel 0.89844 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Si1 core P2  shel P2  shel 0.89844 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Cl2 core P2  shel P2  shel 0.86564 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Si1 core O4  shel P2  shel 0.89844 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Si1 core O4  shel O4  shel 0.89844 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Cl2 core O4  shel P2  shel 0.86564 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Cl2 core O4  shel P3  shel 0.87080 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Cl2 core P2  shel P3  shel 0.87080 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Cl2 core O4  shel O4  shel 0.86564 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Al1 core P2  shel P2  shel 0.87410 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "Al1 core P2  shel P3  shel 0.84520 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "O1  shel Si1 core Si1 core 5.62    163.4  2.4 2.4 5.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "O4  shel Si1 core Cl2 core 5.62    163.4  2.4 2.4 5.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "O4  shel Cl2 core Cl2 core 6.70    166.2  2.4 2.4 5.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "P2  shel Si1 core Al1 core 6.72    162.4  2.4 2.4 5.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "P2  shel Cl2 core Al1 core 6.57    163.4  2.4 2.4 5.0"
    write(10,'(a5)') "three"
    write(10,'(a53)') "P3  shel Cl2 core Al1 core 6.44    150.3  2.4 2.4 5.0"
    write(10,'(a1)') " "
    write(10,'(a10)') "maxcyc 400"
    write(10,'(a1)') " "
    write(10,'(a19)') "gtol opt   9.000000"
    write(10,'(a19)') "ftol opt   9.000000"
    write(10,'(a1)') " "
    write(10,'(a21)') "dump gulp_restart.ff3"
    
    close(10)
                   
  end subroutine gulp_inp_ext_ff3

  subroutine gulp_inp_ssd()

    integer :: ij,ji
    character*3 :: atname
    real(8) :: rcs(3), dcs
    real(8) :: c_s=0.5d0

    open(10,file="gulp_inp.ssd")

    write(10,'(a35)') "opti prop conp phon distance angles"
    write(10,'(a5)') "title"
    write(10,'(a59)') "H-zeolite (correct the title if desirable) Sierka-Sauer DFT"
    write(10,'(a3)') "end"
    write(10,'(a4)') "cell"
    write(10,'(6f11.6)') latt%a
    write(10,'(a10,i5)') "fractional", N_at

    do ij=1,N_at
       if(atoms_fra(ij)%c_s == "c") then
          select case(trim(atoms_fra(ij)%name))
          case ("Al")
             atname="Al1"
          case ("Si")
             atname="Si1"
          case ("F")
             atname="O2"
          case ("O")
             atname="O1"
          case ("H")
             atname="H1"
          end select
          write(10,'(a3,1x,a4,1x,3F12.7)') atname, "core",atoms_fra(ij)%r
       elseif(atoms_fra(ij)%c_s == "s") then
          do ji=1,N_at
             if(atoms_fra(ji)%c_s == "s") cycle
             if(ji==ij) cycle
             rcs=atoms_car(ij)%r-atoms_car(ji)%r
             call image(rcs)
             dcs=sqrt(dot_product(rcs,rcs))
             if(dcs > c_s) cycle
             exit
          end do
          select case(trim(atoms_fra(ji)%name))
          case ("Al")
             atname="Al1"
          case ("Si")
             atname="Si1"
          case ("F")
             atname="O2"
          case ("O")
             atname="O1"
          case ("H")
             atname="H1"
          end select
          write(10,'(a3,1x,a4,1x,3F12.7)') atname, "shel",atoms_fra(ij)%r
       end if
    end do

    write(10,'(a5)') "space"
    write(10,'(a1)') "1"
    write(10,'(a1)') " "
    write(10,'(a7)') "species"
    write(10,'(a13)') "Al1 core  3.0"
    write(10,'(a13)') "Si1 core  4.0"
    write(10,'(a17)') "O1  core  1.22858"
    write(10,'(a17)') "O1  shel -3.22858"
    write(10,'(a17)') "O2  core  0.81753"
    write(10,'(a17)') "O2  shel -2.81753"
    write(10,'(a13)') "H1  core  1.0"
    write(10,'(a1)') " "
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Si1 shel O1 shel 1612.45920    0.29955  0.0      0.0 10.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Si1 shel O2 shel  997.88097    0.33212  0.0      0.0 10.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Al1 shel O1 shel 1395.77463    0.30449  0.0      0.0 10.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "Al1 shel O2 shel 1644.88177    0.29139  0.0      0.0 10.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "H1  shel O1 shel 7614.58003    0.19913  0.0      0.0 10.0"
    write(10,'(a4)') "buck"
    write(10,'(a57)') "H1  shel O2 shel  368.64803    0.22511  0.0      0.0 10.0"
    write(10,'(a1)') " "
    write(10,'(a6)') "spring"
    write(10,'(a12)') "O1 122.47853"
    write(10,'(a12)') "O2  70.15123"
    write(10,'(a1)') " "
    write(10,'(a5)') "three"
    write(10,'(a52)') "Si1 core O1 shel O1 shel 0.144703 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a52)') "Si1 core O1 shel O2 shel 0.384711 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a52)') "Al1 core O1 shel O1 shel 0.893930 109.47 2.4 2.4 4.0"
    write(10,'(a5)') "three"
    write(10,'(a52)') "Al1 core O1 shel O2 shel 0.686678 109.47 2.4 2.4 4.0"
    write(10,'(a1)') " "
    write(10,'(a10)') "maxcyc 400"
    write(10,'(a1)') " "
    write(10,'(a19)') "gtol opt   9.000000"
    write(10,'(a19)') "ftol opt   9.000000"
    write(10,'(a1)') " "
    write(10,'(a21)') "dump gulp_restart.ssd"
    
    close(10)

  end subroutine gulp_inp_ssd

end program make_zeo

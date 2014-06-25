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
program genlat

  !This program serves to generate Na*Nb*Nc super cell
  !from any lattice unit cell.
  !Input data can be specified in both cartesian or fractional
  !formats
  !Output coordinates also can be obtained in both formats
  !Program can be used as the converter coordinates
  !
  !Example of input file:
  !cell (vect)
  !a b c alpha beta gamma       - unit cell parameters
  !or
  !x1 y1 z1
  !x2 y2 z2
  !x3 y3 z3
  !frac (cart)                  - type of input coordinates
  !N                            - number of atoms into the unit cell
  !Atom1    x1     y1     z1 \
  !.........................  | - atom names and coordinates
  !AtomN    xN     yN     zN /
  !Na Nb Nc                     - translation along a, b and c directions
  !xd yd zd                     - displasement along a, b and c directions
  !                               in frac or cart (depends on coordinate type)
  !frac (cart)                  - type of output coordinates
  !sym
  !
  !Execution : genlat < input_file > output_file
  !
  implicit none

  integer*4 :: natoms,na,nb,nc
  real*8 :: xd,yd,zd
  integer*4 :: i,j,k,l,n
  real*8 :: cart(3,10000),frac(3,10000) !,fbuf(3,10000)
  character*4 :: atnm(10000)
  real*8 :: bmat(3,3),v1(3),v2(3),v3(3),x1(3),y1(3),z1(3)
  real*8 :: a,b,c,alpha,beta,gamma,a1,b1,c1
  character*4 :: coor_type,save_type
  character*4 :: type_lat
  character*3 :: sym
  logical :: do_sym
  real*8, parameter :: pi=3.141592653589793d0, pi_degree=180.0d0, a2b=1.889725988579d0

  do_sym=.false.
  frac=0.0d0; cart=0.0d0

  read(5,*) type_lat
  if(type_lat == "cell") then
       read(5,*) a,b,c,alpha,beta,gamma
  elseif(type_lat == "vect") then
       read(5,*) v1
       read(5,*) v2
       read(5,*) v3
  else
       stop "wrong keyword"
  endif
  read(5,*) coor_type
  if(coor_type /= "frac" .and. coor_type /= "cart") &
       stop 'wrong coordinate type'
  read(5,*) natoms
  do i=1,natoms
     if(coor_type=="frac") then
        read(5,*) atnm(i), frac(:,i)
        do j=1,3
           if(frac(j,i) < 0.0d0) then 
              frac(j,i)=1.0d0+frac(j,i)
           elseif(frac(j,i) > 1.0d0) then
              frac(j,i)=frac(j,i)-1.0d0
           end if
        end do
     else if(coor_type=="cart") then
        read(5,*) atnm(i), cart(:,i)
     end if
  end do
  read(5,*) na,nb,nc
  read(5,*) xd,yd,zd
  read(5,*) save_type
  read(5,*,end=100) sym
  if(sym=="sym") do_sym=.true.
100 continue 

  if(type_lat == "cell") call cell2vect(a,b,c,alpha,beta,gamma,v1,v2,v3)
  if(type_lat == "vect") call vect2cell(v1,v2,v3,a,b,c,alpha,beta,gamma)

  if(coor_type=="cart") call cart2frac()

!  fbuf(1:3,1:natoms)=frac(1:3,1:natoms)
  n=0
  l_na: do i=0,na-1
     l_nb: do j=0,nb-1
        l_nc: do k=0,nc-1
           l_ns: do l=1,natoms
              n=n+1
              frac(1,n)=frac(1,l)+i
              frac(2,n)=frac(2,l)+j
              frac(3,n)=frac(3,l)+k
              atnm(n)=atnm(l)
           end do l_ns
        end do l_nc
     end do l_nb
  end do l_na

  !converting fractional coordinates to cartesian
  natoms=n
!  do l=1,natoms
!     frac(1,l)=frac(1,l)+xd
!     frac(2,l)=frac(2,l)+yd
!     frac(3,l)=frac(3,l)+zd
!  enddo 
  call frac2cart()

  a1=a; b1=b; c1=c
  a=a*na; b=b*nb; c=c*nc

  call cell2vect(a,b,c,alpha,beta,gamma,v1,v2,v3)
  call cart2frac()

  if(coor_type=="cart") then
     do l=1,natoms
        cart(1,l)=cart(1,l)+xd
        cart(2,l)=cart(2,l)+yd
        cart(3,l)=cart(3,l)+zd
     enddo
  elseif(coor_type=="frac") then
     call fd2cd(xd,yd,zd)
     do l=1,natoms
        cart(1,l)=cart(1,l)+xd
        cart(2,l)=cart(2,l)+yd
        cart(3,l)=cart(3,l)+zd
     enddo
  endif

  if(do_sym .and. (abs(xd)<0.001d0 .and. abs(yd)<0.001d0 .and. abs(zd)<0.001d0)) then
     do l=1,natoms
        if(na > 1) frac(1,l)=frac(1,l)-0.5d0
        if(nb > 1) frac(2,l)=frac(2,l)-0.5d0
        if(nc > 1) frac(3,l)=frac(3,l)-0.5d0
     enddo
     call frac2cart()
  endif

  write(6,'(a7)')  " $title"
  write(6,'(a11)') " !!!!!!!!!!"
  if(save_type=='frac') then 
     write(6,'(a20)')  " $coord fractional 1"
  else if(save_type=='cart') then
     write(6,'(a9)')   " $coord 1"
  end if
  do i=1,natoms
     if(save_type=='cart') then
        write(6,'(3f16.9,1x,a6)')  cart(:,i),atnm(i)
     else if(save_type=='frac') then
        write(6,'(3f16.9,1x,a6)') frac(:,i),atnm(i)
     end if
  end do
  write(6,'(a11,6f14.8)') " $unitcell ",a,b,c,alpha,beta,gamma
  write(6,'(a5)')  " $end"
  write(6,'(a1)') " "
  write(6,'(3f14.8,a3,3f14.8)') v1,' | ',v1*a2b
  write(6,'(3f14.8,a3,3f14.8)') v2,' | ',v2*a2b
  write(6,'(3f14.8,a3,3f14.8)') v3,' | ',v3*a2b
  write(6,*) "*************************************************"
  write(6,*) "*************************************************"
  write(6,'(i7)') natoms
  write(6,'(a7)') "Lattice"
  do i=1,natoms
     write(6,'(a6,1x,3f16.9,a3,3f16.9)')  atnm(i),cart(:,i),' | ',cart(:,i)*a2b
  enddo
contains
  !******************************************************************
  subroutine cell2vect(ac,bc,cc,alphac,betac,gammac,vc1,vc2,vc3)

    real*8, intent(in) :: ac,bc,cc,alphac,betac,gammac
    real*8, intent(out) :: vc1(3),vc2(3),vc3(3)

    real*8 :: alpha_r, beta_r, gamma_r

    alpha_r=pi*alphac/pi_degree
    beta_r=pi*betac/pi_degree
    gamma_r=pi*gammac/pi_degree

    vc1(1)=ac; vc1(2)=0.0d0; vc1(3)=0.0d0
    vc2(1)=bc*cos(gamma_r); vc2(2)=bc*sin(gamma_r); vc2(3)=0.0d0
    vc3(1)=cc*cos(beta_r)
    vc3(2)=(bc*cc*cos(alpha_r)-vc2(1)*vc3(1))/vc2(2)
    vc3(3)=sqrt(cc*cc-vc3(1)*vc3(1)-vc3(2)*vc3(2))

  end subroutine cell2vect
  !******************************************************************
  subroutine vect2cell(vc1,vc2,vc3,ac,bc,cc,alphac,betac,gammac)

    real*8, intent(in) :: vc1(3),vc2(3),vc3(3)
    real*8, intent(out) :: ac,bc,cc,alphac,betac,gammac
    real*8 :: alpha_r, beta_r, gamma_r
    integer*4 :: i

     ac=sqrt(dot_product(vc1,vc1))
     bc=sqrt(dot_product(vc2,vc2))
     cc=sqrt(dot_product(vc3,vc3))

     alpha_r=0.0d0; beta_r=0.0d0; gamma_r=0.0d0;
     do i=1,3
        alpha_r=alpha_r+vc2(i)*vc3(i)
        beta_r=beta_r+vc1(i)*vc3(i)
        gamma_r=gamma_r+vc2(i)*vc1(i)
     end do

     alpha_r=acos(alpha_r/(bc*cc))
     beta_r=acos(beta_r/(ac*cc))
     gamma_r=acos(gamma_r/(ac*bc))

     alphac=pi_degree*alpha_r/pi
     betac=pi_degree*beta_r/pi
     gammac=pi_degree*gamma_r/pi

  end subroutine vect2cell
 !******************************************************************
  subroutine cart2frac()

    real*8 :: bmata(3,3),det,rr
    integer*4 :: i,j

    bmata(1,1)=v2(2)*v3(3)-v2(3)*v3(2)
    bmata(2,1)=v1(3)*v3(2)-v1(2)*v3(3)
    bmata(3,1)=v1(2)*v2(3)-v1(3)*v2(2)
    bmata(1,2)=v2(3)*v3(1)-v2(1)*v3(3)
    bmata(2,2)=v1(1)*v3(3)-v1(3)*v3(1)
    bmata(3,2)=v1(3)*v2(1)-v1(1)*v2(3)
    bmata(1,3)=v2(1)*v3(2)-v2(2)*v3(1)
    bmata(2,3)=v1(2)*v3(1)-v1(1)*v3(2)
    bmata(3,3)=v1(1)*v2(2)-v1(2)*v2(1)

    det=v1(1)*bmata(1,1)+v2(1)*bmata(2,1)+v3(1)*bmata(3,1)
    rr=0.0d0
    if(abs(det) > 0.0d0) rr=1.0d0/det

    bmata=bmata*rr

    do i=1,natoms
       frac(:,i)=matmul(bmata,cart(:,i))
       do j=1,3
          if(abs(frac(j,i)) /= 0.5d0) then
             frac(j,i)=frac(j,i)-anint(frac(j,i))
          end if
          if(frac(j,i) < 0.0d0) then 
             frac(j,i)=1.0d0+frac(j,i)
          elseif(frac(j,i) > 1.0d0) then
             frac(j,i)=frac(j,i)-1.0d0
          end if
       end do
    end do

  end subroutine cart2frac
  !******************************************************************
  subroutine frac2cart()

    real*8 :: bmata(3,3)
    integer*4 :: i

    bmata(:,1)=v1
    bmata(:,2)=v2
    bmata(:,3)=v3

    do i=1,natoms
       cart(:,i)=matmul(bmata,frac(:,i))
    end do

  end subroutine frac2cart
  !******************************************************************
  subroutine fd2cd(x,y,z)

    real*8 :: x,y,z
    real*8 :: bmata(3,3),r(3),r1(3)
    integer*4 :: i

    bmata(:,1)=v1
    bmata(:,2)=v2
    bmata(:,3)=v3
 
    r(1)=x; r(2)=y; r(3)=z

    r1(:)=matmul(bmata,r(:))

    x=r1(1); y=r1(2); z=r1(3)

  end subroutine fd2cd
  !******************************************************************
end program genlat

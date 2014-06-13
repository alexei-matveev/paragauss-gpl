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
module mesh_module

 implicit none
 
 type sphere_mesh
 real(kind=8), dimension(110) :: xp
 real(kind=8), dimension(110) :: yp
 real(kind=8), dimension(110) :: zp
 real(kind=8), dimension(110) :: w
 integer(kind=4)              :: npoints  !  number of points on sphere
 end type sphere_mesh
 
 type(sphere_mesh)               mesh
 real(kind=8), parameter      :: pi=3.14159265358979_8        

contains 

   subroutine genmesh(rc,rsphere,ispher)
   
   implicit none
   integer(kind=4)               npoints,ispher, l, m
   real   (kind=8)               a1(6,3),  a2(12,3),    &
                                 a3(8,3),  b1(24,3),    &
                                 b2(24,3), b3(24,3),    & 
                                 c1(24,3), k(3,3),      &
                                 mk(3,3),  pk(3),qk(3), &
                                 lk(3,3),  rc(3),       & ! rc - centre of sphere 
                                 rsphere,               & ! radius 
                                 r1, r2,                &                                          
                                 asf(110,3)               ! weights fixed                      
   integer(kind=4),dimension(3)  :: nspnt                 ! 3 schemes are available

   intent(in)  ::  rc, rsphere, ispher
                                                              
 data lk /   0.00000000000000d0,    0.00000000000000d0,   0.00000000000000d0, &
             0.30151134457800d0,    0.00000000000000d0,   0.00000000000000d0, &
             0.18511563534500d0,    0.39568947305600d0,   0.69042104838200d0  /
                                                                                 
 data mk /   0.00000000000000d0,    0.00000000000000d0,   0.00000000000000d0, &
             0.90453403373300d0,    0.00000000000000d0,   0.00000000000000d0, &
             0.96512403508700d0,    0.82876998125300d0,   0.21595729184600d0  /
 
 data pk /   0.88807383300000d0,    0.00000000000000d0,   0.87815891060400d0  /

 data qk /   0.45970084300000d0,    0.00000000000000d0,   0.47836902881200d0  /

 data asf/ 6*0.00952380952400d0,  8*0.03214285714300d0,24*0.02857142857100d0, &
          72*0.00000000000000d0,  6*0.01269841269800d0,12*0.02257495590800d0, &
           8*0.02109375000000d0, 24*0.02017333553800d0,60*0.00000000000000d0, &
           6*0.00382827049494d0,  8*0.00979373751249d0,24*0.00821173728319d0, &
          24*0.00959547133607d0, 24*0.00994281489118d0,24*0.00969499636166d0  /
                                                                              
 data nspnt / 38, 50, 110 /    

! = = = = = = = = = = = = = a1 = = = = = = = = = = = = = = =

      l=1
      call sfc(0.d0,0.d0,1.d0,0,0,1,a1,6,l)
      call sfc(0.d0,1.d0,0.d0,0,1,0,a1,6,l)
      call sfc(1.d0,0.d0,0.d0,1,0,0,a1,6,l)

! = = = = = = = = = = = = = a2 = = = = = = = = = = = = = = = =

      l=1
      r1=1.d0/sqrt(2.d0)
      call sfc(r1,r1,0.d0,1,1,0,a2,12,l)
      call sfc(r1,0.d0,r1,1,0,1,a2,12,l)
      call sfc(0.d0,r1,r1,0,1,1,a2,12,l)

! = = = = = = = = = = = = = a3 = = = = = = = = = = = = = = = =

      l=1
      r1=1.d0/sqrt(3.d0)
      call sfc(r1,r1,r1,1,1,1,a3,8,l)

! = = = = = = = = = = = = = b1 = = = = = = = = = = = = = = = =

      l=1
      r1=lk(1,ispher)
      r2=mk(1,ispher)
      call sfc(r1,r1,r2,1,1,1,b1,24,l)
      call sfc(r1,r2,r1,1,1,1,b1,24,l)
      call sfc(r2,r1,r1,1,1,1,b1,24,l)

! = = = = = = = = = = = = = b2 = = = = = = = = = = = = = = = =

      l=1
      r1=lk(2,ispher)
      r2=mk(2,ispher)
      call sfc(r1,r1,r2,1,1,1,b2,24,l)
      call sfc(r1,r2,r1,1,1,1,b2,24,l)
      call sfc(r2,r1,r1,1,1,1,b2,24,l)

! = = = = = = = = = = = = = b3 = = = = = = = = = = = = = = = =

      l=1
      r1=lk(3,ispher)
      r2=mk(3,ispher)
      call sfc(r1,r1,r2,1,1,1,b3,24,l)
      call sfc(r1,r2,r1,1,1,1,b3,24,l)
      call sfc(r2,r1,r1,1,1,1,b3,24,l)

! = = = = = = = = = = = = = c1 = = = = = = = = = = = = = = = =

      l=1
      r1=pk(ispher)
      r2=qk(ispher)
      call sfc(r1,r2,0.d0,1,1,0,c1,24,l)
      call sfc(r1,0.d0,r2,1,0,1,c1,24,l)
      call sfc(0.d0,r1,r2,0,1,1,c1,24,l)
      call sfc(r2,r1,0.d0,1,1,0,c1,24,l)
      call sfc(r2,0.d0,r1,1,0,1,c1,24,l)
      call sfc(0.d0,r2,r1,0,1,1,c1,24,l)

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 
      npoints = 0
      do  m = 1, nspnt(ispher)
        mesh%w(m) = 4*pi*asf(m,ispher)
      end do

                call coor(rsphere,rc,a1, 6)     
      select case(ispher)
        case(1)
                call coor(rsphere,rc,a3, 8)
                call coor(rsphere,rc,c1,24)     
        case(2)
                call coor(rsphere,rc,a2,12)
                call coor(rsphere,rc,a3, 8)
                call coor(rsphere,rc,b1,24)
        case(3)
                call coor(rsphere,rc,a3, 8)
                call coor(rsphere,rc,b1,24)
                call coor(rsphere,rc,b2,24) 
                call coor(rsphere,rc,b3,24) 
                call coor(rsphere,rc,c1,24) 
        end select
                               
     if(mesh%npoints .ne. nspnt(ispher) ) then
        write(*,*) 'The number of generated and prescribed points is not equal !!!'
        stop
     end if
     
     write(*,*) 'out mesh'
     return
   end subroutine genmesh
  
   subroutine sfc(ax,ay,az,ix,iy,iz,a,ndim,l)
   
   implicit none
   integer(kind=4), intent(in)        :: ix,iy,iz,ndim
   integer(kind=4), intent(inout)     :: l
   integer(kind=4)                       isgnx,isgny,isgnz,nx,ny,nz
   integer(kind=4)                       ipx,ipy,ipz
   real(kind=8),    dimension(ndim,3) :: a
   real(kind=8),    intent(in)        :: ax,ay,az
      isgnx=1
      isgny=1
      isgnz=1
      nx=2
      ny=2
      nz=2
      if(ix.eq.0) nx=1
      if(iy.eq.0) ny=1
      if(iz.eq.0) nz=1
      do ipx=1,nx
         do ipy=1,ny
            do ipz=1,nz
            a(l,1)=ax*isgnx
            a(l,2)=ay*isgny
            a(l,3)=az*isgnz
            if(iz.ne.0) isgnz=-isgnz
            l=l+1
            end do
         if(iy.ne.0) isgny=-isgny
         end do
      if(ix.ne.0) isgnx=-isgnx
      end do
      return
      end subroutine sfc 

   subroutine coor(rad,rc,a,m)

   implicit none
   integer(kind=4)                  i,m 
   real(kind=8), dimension(m,3)  :: a
   real(kind=8), dimension(3)    :: rc
   real(kind=8)                     rad

      do i=1,m
      mesh%npoints=mesh%npoints+1
      mesh%xp(mesh%npoints) = rc(1) + rad*a(i,1)
      mesh%yp(mesh%npoints) = rc(2) + rad*a(i,2)
      mesh%zp(mesh%npoints) = rc(3) + rad*a(i,3)
      end do
      return
   end subroutine coor

end module mesh_module
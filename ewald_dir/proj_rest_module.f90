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
module proj_rest_module
use mesh_module
 
  implicit none
  real(kind=8),              dimension(400)  :: vlm       !   projected potential
  integer(kind=4)                               nlm       !   No components

contains

subroutine sh(x,y,z,lmax,c)
!
!   Recursive calculations for Normalized Real Solid Harmonics / r**l
!
   implicit none
   
   real(kind=8), dimension(:,:) :: c ! ylm(1:lmax+1, 1:2*lmax+1)
   real(kind=8)                    x,y,z,r,sqr
   integer(kind=4)                 lmax,l,m,im,i
   
   c(1,1) = 1.d0 
   if(lmax==0) return
   r = sqrt(x**2 +y**2 + z**2)
   
    if(r==0.d0) then
     c(2:lmax+1,:)=0.d0
     return
    end if
    
   c(2,1) = z / r
   c(2,2) = x / r
   c(2,3) = y / r
   if(lmax==1) return 
 
   do l=3,lmax+1
   m=1
   c(l,m)=p(l-1,z/r)
 
   do im=1,l-3
      do i=1,2
        m=m+1   
        c(l,m) = 1.d0/sqrt( float( (l-1)**2-(m/2)**2) ) *            &
              ( (2*l-3) * z/r * c(l-1,m) -                           &
              sqrt( float( (l-2)**2-(m/2)**2 ) ) * c(l-2,m) )  
      end do 
   end do  ! im
        
   im=l-2
   do i=1,2
     m=m+1
     c(l,m) = 1.d0/sqrt( float( (l-1)**2-(m/2)**2) ) *                & 
                (2*l-3) * z/r * c(l-1,m)         
   end do
        
   im=l-1
   sqr=sqrt( float(2*l-3)/float(2*l-2) )
   m=m+1
   c(l,m) = sqr * ( x*c(l-1,2*l-4) - y*c(l-1,2*l-3) ) / r                       
   m=m+1
   c(l,m) = sqr * ( y*c(l-1,2*l-4) + x*c(l-1,2*l-3) ) / r 
                       
   end do  ! l
   return

contains
  
  function p(n,t)
  implicit none
  real(kind=8)    p,t,pk,pk1,pk2
  integer(kind=4) n,k

  select case(n)
    case(0)
      p=1.d0
      return
    case(1)
      p=t
      return
  end select
      pk2=1.d0
      pk1=t
        do  k=2,n
          pk=( (2*k-1)*t*pk1 - (k-1)*pk2 ) / k
          pk2=pk1
          pk1=pk
        end do
      p=pk
      return
  end function p

end subroutine sh  
!
!        lmax=0,1,2,...for (s,p,d,...) contributions         
!  it would be more convenient to rearrange loop 2 and 1 ... 
!
   subroutine project(vs,lmax,rc,c)      
   
   implicit none
   real(kind=8), dimension(:,:) :: c 
   real(kind=8)                    vs(:)   ! Input potential on sphere
   real(kind=8)                    rc(3)
   integer(kind=4) l,lmax,m,ip
   intent(in)   :: lmax 
 
    vlm(:)=0.d0
    do ip=1,mesh%npoints
      call sh(mesh%xp(ip)-rc(1),mesh%yp(ip)-rc(2),mesh%zp(ip)-rc(3),lmax,c)
      nlm=0
      do l=1,lmax+1
        do m=1,2*l-1
        nlm=nlm + 1
           vlm(nlm) = vlm(nlm) + vs(ip) * c(l,m) * mesh%w(ip)        
        end do
      end do
   end do

   write(6,'(/30(1h*), " Y_lm amlitudes ",20(1h*)/t10,"N", t15,"l", t20,"m", t35,"Vlm")' )
   nlm=0
   do l=1,lmax+1
     do m=1,2*l-1
     nlm=nlm+1
      write(6,'( t8,i3,t14,i2,t19,i2,t27,f19.15,3x)') & 
      nlm,l,m,vlm(nlm)
     end do
   end do   
   write(*,'(/67(1h*)/)')
      
   return
   end subroutine project
!
!     Restore potential in (x,y,z)-point inside a sphere
!
   subroutine rest(x,y,z,rs,rc,lmax,c,vrest)

   implicit none
   real(kind=8), dimension(:,:) :: c 
   real(kind=8), dimension(3)   :: rc
   real(kind=8)                    x,y,z,rs,r,rrl,vrest
   integer(kind=4)                 lm,lmax,l,m
   intent(in)                      x,y,z,rs,lmax


      call sh(x-rc(1),y-rc(2),z-rc(3),lmax,c)    
      vrest = 0.d0
      rrl=1.d0
      r=sqrt( (x-rc(1))**2 + (y-rc(2))**2 + (z-rc(3))**2 )
      lm=0
        do l=1,lmax+1
        if(l .gt. 1) rrl = rrl * (r/rs) 
          do m=1,2*l-1
          lm=lm+1
          vrest = vrest + rrl * float(2*l-1)/(4*pi) * c(l,m)*vlm(lm)
          end do
        end do
      return
   end subroutine rest 

end module proj_rest_module

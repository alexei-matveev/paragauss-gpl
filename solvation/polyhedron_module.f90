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
!===============================================================
! Public interface of module
!===============================================================
module polyhedron_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !  Author: AS
  !  Date: 03/2010
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

# include <def.h>
  use type_module ! type specification parameters
  use datatype
  use iounitadmin_module, only: output_unit
  use output_module ,only: output_cavity_data
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------
  type, public :: triangles
     real(r8_kind),pointer    :: xyz(:,:)
     real(r8_kind),pointer    :: xyz_centers(:,:)
     integer(i4_kind),pointer :: index(:,:)
     real(r8_kind)            :: radius
     real(r8_kind)            :: area
  end type triangles

  !------------ Declaration of constants and variables ------------
  type(triangles) ,public,target    :: surf_elem    
  integer(i4_kind),public           :: N_points_of_triangles
  integer(i4_kind),public           :: N_centers_on_sphere
  integer(i4_kind),parameter,public :: MAX_ATOMS_QMMM=40
  integer(i4_kind),parameter,public :: MAX_ATOMS_QMMM1=100
  integer(i4_kind),public           :: NDIV=5
  integer(i4_kind),public           :: n_rotations,max_cent
  real(r8_kind), public             :: radius
  integer(i4_kind), public          :: local_point_factor

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public generate_octahedron,generate_cube,generate_dodecahedron, &
         generate_doublepyramide,more_triangles

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------
  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine generate_octahedron(N_atom_spheres,do_cavitation,do_gradients)
    ! Purpose: generate set of surface points on a sphere starting from a cube
    !          and do subdivisions
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),intent(in) :: N_atom_spheres
    logical,intent(in)          :: do_cavitation, do_gradients
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind),pointer     :: xyz(:,:),xyz_cent(:,:)
    integer(i4_kind), pointer :: indexx(:,:)
    real(r8_kind)             :: a,d
    real(r8_kind)             :: xyz_buf(3)
    integer(i4_kind)          :: N_dim_xyz,N_dim_cent,N_dim_xyz_st,N_dim_cent_st
    integer(i4_kind)          :: N_dim_xyz_next,N_dim_cent_next
    integer(i4_kind)          :: i, j, l1, n1, i1, status1
    !------------ Executable code --------------------------------
    
    a=1.0_r8_kind
    radius=a
    d=sqrt(2.0_r8_kind)*a
    
    N_dim_cent=8
    N_dim_cent_st=8
    N_dim_xyz=6
    N_dim_xyz_st=6

    local_point_factor=2
    if(N_atom_spheres >= MAX_ATOMS_QMMM1) local_point_factor=1
    
    i=1
    do
       if(i>local_point_factor-1) exit
       N_dim_xyz=N_dim_xyz+(N_dim_cent*3)/2
       N_dim_cent=N_dim_cent*4
       if(i==1) then
          N_dim_xyz_next=N_dim_xyz
          N_dim_cent_next=N_dim_cent
       endif
       i=i+1
    enddo
    
    allocate(surf_elem%xyz(N_dim_xyz,3),stat=status1)
    ASSERT(status1==0)
    
    !definition of initial set of points
    xyz=>surf_elem%xyz
    xyz=0.0_r8_kind

    xyz(1,1)=a
    xyz(2,1)=-a
    xyz(3,2)=a
    xyz(4,2)=-a
    xyz(5,3)=a
    xyz(6,3)=-a
    
    !definition of triangles

    allocate(surf_elem%index(N_dim_cent,3),stat=status1)
    ASSERT(status1==0)

    indexx=>surf_elem%index

    indexx(1,1)=1
    indexx(1,2)=3
    indexx(1,3)=5
    
    indexx(2,1)=1
    indexx(2,2)=4
    indexx(2,3)=5

    indexx(3,1)=1
    indexx(3,2)=3
    indexx(3,3)=6

    indexx(4,1)=1
    indexx(4,2)=4
    indexx(4,3)=6

    indexx(5,1)=2
    indexx(5,2)=3
    indexx(5,3)=5

    indexx(6,1)=2
    indexx(6,2)=4
    indexx(6,3)=5

    indexx(7,1)=2
    indexx(7,2)=3
    indexx(7,3)=6
    
    indexx(8,1)=2
    indexx(8,2)=4
    indexx(8,3)=6

    do j=1,local_point_factor-1
#ifdef _LINUX1
       call more_triangles(N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#else
       call more_triangles(xyz,indexx,N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#endif
       N_dim_xyz_st=N_dim_xyz_next
       N_dim_cent_st=N_dim_cent_next
       N_dim_xyz_next= N_dim_xyz_st+(N_dim_cent_st*3)/2
       N_dim_cent_next= N_dim_cent_st*4
    enddo

    ! definition of triangle centers
    allocate(surf_elem%xyz_centers(N_dim_cent,3),stat=status1)
    ASSERT(status1==0)

    xyz_cent=>surf_elem%xyz_centers
    do i=1,N_dim_cent_st
       i1=indexx(i,1)
       l1=indexx(i,2)
       n1=indexx(i,3)
       xyz_buf=(xyz(i1,:)+xyz(l1,:)+xyz(n1,:))/3.0_r8_kind
       xyz_cent(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
    enddo
    N_points_of_triangles=N_dim_xyz_st
    N_centers_on_sphere=N_dim_cent_st
    
    if (do_cavitation .and. output_cavity_data .and. .not.do_gradients) &
         write (output_unit,'(a42,i5)') &
         'the number of triangles on each sphere is ', N_centers_on_sphere
    
  end subroutine generate_octahedron
  !*************************************************************

  !*************************************************************
  subroutine generate_cube(point_factor,do_cavitation,do_gradients)
    ! Purpose: generate set of surface points on a sphere starting from a cube
    !          and do subdivisions
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),intent(in) :: point_factor
    logical,intent(in)          :: do_cavitation,do_gradients
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind),pointer     :: xyz(:,:),xyz_cent(:,:)
    integer(i4_kind), pointer :: indexx(:,:)
    real(r8_kind)             :: a,d
    real(r8_kind)             :: r_buf,xyz_buf(3)
    integer(i4_kind)          :: tang_buf(4)
    integer(i4_kind)          :: N_dim_xyz,N_dim_cent,N_dim_xyz_st,N_dim_cent_st
!!!MF >>>>
    integer(i4_kind)          :: N_dim_xyz_next,N_dim_cent_next
!!!MF <<<<
    integer(i4_kind)          :: i,j,k,l,m,n,l1,n1,ind,i1,status1
    !------------ Executable code --------------------------------

    a=sqrt(2.0_r8_kind)/2.0_r8_kind
    radius=sqrt(1.5_r8_kind)
    d=(radius-a)**2+1.0_r8_kind

    N_dim_cent=24
    N_dim_cent_st=24
    N_dim_xyz=14
    N_dim_xyz_st=14

    if (do_cavitation) then
       local_point_factor=0
    else
       local_point_factor=point_factor
    endif

!!!MF >>>>
    i=1
    do 
       if(i>local_point_factor .and. & 
            (N_dim_cent>=max_cent  .or. do_cavitation) ) exit
       N_dim_xyz=N_dim_xyz+N_dim_cent*3/2
       N_dim_cent=N_dim_cent*4
       if(i==1) then
          N_dim_xyz_next=N_dim_xyz
          N_dim_cent_next=N_dim_cent
       endif
       i=i+1
    enddo
!!!MF <<<<
!     if(local_point_factor < i) local_point_factor=i-1
!     if(do_cavitation) local_point_factor=0

    allocate(surf_elem%xyz(N_dim_xyz,3),stat=status1)
    ASSERT(status1 == 0)

    !definition of initial set of points
    xyz=>surf_elem%xyz
    xyz=0.0_r8_kind

    xyz(1:4,1)=a
    xyz(5:8,1)=-a
    xyz(9,1)=radius
    xyz(10,1)=-radius
    xyz(11:14,1)=0.0_r8_kind

    xyz(1,2)=-a
    xyz(2,2)=a
    xyz(3,2)=-a
    xyz(4,2)=a
    xyz(5,2)=-a
    xyz(6,2)=a
    xyz(7,2)=-a
    xyz(8,2)=a
    xyz(9:10,2)=0.0_r8_kind
    xyz(11,2)=radius
    xyz(12,2)=-radius
    xyz(13:14,2)=0.0_r8_kind

    xyz(1:2,3)=a
    xyz(3:4,3)=-a
    xyz(5:6,3)=a
    xyz(7:8,3)=-a
    xyz(9:12,3)=0.0_r8_kind
    xyz(13,3)=radius
    xyz(14,3)=-radius

    !definition of triangles

    allocate(surf_elem%index(N_dim_cent,3),stat=status1)
    ASSERT(status1 == 0)

    surf_elem%index=0

    indexx=>surf_elem%index

    ind=1
    do i=1,6
       k=i+8
       m=4*(i-1)+1
       indexx(m:m+3,1)=k
       
       n=1
       do j=1,8
          xyz_buf=xyz(k,:)-xyz(j,:)
          r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
          if (r_buf <= sqrt(d)+0.01_r8_kind) then
             tang_buf(n)=j
             n=n+1
             if (n > 4) exit
          endif

       enddo

       do l=1,3
          l1=tang_buf(l)
          do n=l+1,4
             n1=tang_buf(n)
             xyz_buf=xyz(n1,:)-xyz(l1,:)
             r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
             if (r_buf <= 2*a+0.01_r8_kind) then
                indexx(ind,2)=l1
                indexx(ind,3)=n1
                ind=ind+1
             endif
          enddo
       enddo
    enddo

    !dividing each triangles in four new triangles 
!!!MF >>>>
    do j=1,local_point_factor
#ifdef _LINUX1
       call more_triangles(N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#else
       call more_triangles(xyz,indexx,N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#endif
       N_dim_xyz_st=N_dim_xyz_next
       N_dim_cent_st=N_dim_cent_next
       N_dim_xyz_next= N_dim_xyz_st+N_dim_cent_st*3/2
       N_dim_cent_next= N_dim_cent_st*4
    enddo
!!!MF <<<<

DPRINT 'maximum allowd centers',N_dim_cent
    ! definition of triangle centers
    allocate(surf_elem%xyz_centers(N_dim_cent,3),stat=status1)
    ASSERT(status1 == 0)

    xyz_cent=>surf_elem%xyz_centers
    do i=1,N_dim_cent_st
       i1=indexx(i,1)
       l1=indexx(i,2)
       n1=indexx(i,3)
       xyz_buf=(xyz(i1,:)+xyz(l1,:)+xyz(n1,:))/3.0_r8_kind
       xyz_cent(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
    enddo

    N_points_of_triangles=N_dim_xyz_st
    N_centers_on_sphere=N_dim_cent_st

    if (do_cavitation .and. output_cavity_data .and. .not.do_gradients) &
         write (output_unit,'(a42,i5)') 'the number of triangles on each sphere is', N_centers_on_sphere
!     nullify(xyz,xyz_cent,index)
  end subroutine generate_cube
  !*************************************************************

  !*************************************************************
  subroutine generate_dodecahedron(gp93,point_factor,name_point_group, &
       do_cavitation,do_gradients)
    !  Purpose: generate points due to a pentakis-dodecahedron
    !           rotate in right symmetry
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),intent(in) :: gp93
    integer(i4_kind),intent(in) :: point_factor
    logical,intent(in)          :: do_cavitation,do_gradients
    character(len=4),intent(in)      :: name_point_group
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind),pointer     :: xyz(:,:),xyz_cent(:,:)
    integer(i4_kind), pointer :: indexx(:,:)
    real(r8_kind)             :: x,y,z
    real(r8_kind)             :: t,s,h,r,b,c1,c2,h1,h2,rr,cosT,sinT
    real(r8_kind)             :: r_buf,xyz_buf(3)
    integer(i4_kind)          :: tang_buf(5)
    integer(i4_kind)          :: N_dim_xyz,N_dim_cent,N_dim_xyz_st,N_dim_cent_st
!!!MF >>>>
    integer(i4_kind)          :: N_dim_xyz_next,N_dim_cent_next
!!!MF <<<<
!!!MF 12.4.2000 little rotations for D3 and D5 cos/sin of 2 pi/12 or 2 pi/20 >>>>
    real(r8_kind), parameter  :: cos12=0.8660254037844_r8_kind
    real(r8_kind), parameter  :: sin12=0.5_r8_kind
    real(r8_kind), parameter  :: cos20=0.9510565162952_r8_kind
    real(r8_kind), parameter  :: sin20=0.3090169943749_r8_kind
!!!MF <<<<
    integer(i4_kind)          :: i,j,k,l,m,n,l1,n1,ind,i1,status2
    !------------ Executable code --------------------------------

    !definition a pentakis dodecahedron 
    t=(sqrt(5.0_r8_kind)+1.0_r8_kind)/2.0_r8_kind
    s=sqrt(3.0_r8_kind-t)/2.0_r8_kind      
    h=(2.0_r8_kind*s+(t+1.0_r8_kind)*sqrt(t+2.0_r8_kind))/4.0_r8_kind
    radius=sqrt(h**2+s**2)
    b=radius*sqrt(3.0_r8_kind)/3.0_r8_kind
    c1=(2.0_r8_kind*b+h)/5.0_r8_kind
    c2=(2.0_r8_kind*h+2.0_r8_kind*b+s)/5.0_r8_kind
    rr=sqrt(c1**2+c2**2)
    h1=c1*radius/rr
    h2=c2*radius/rr
    cosT=h2/sqrt(h1**2+h2**2)
    sinT=sqrt(1.0_r8_kind-cosT**2)
    r=sqrt(radius**2-1.0_r8_kind)

    N_dim_cent=60
    N_dim_cent_st=60
    N_dim_xyz=32
    N_dim_xyz_st=32

    if (do_cavitation) then
       local_point_factor=1
    else
       local_point_factor=point_factor
    endif

    if(gp93==1) then 
       local_point_factor=NDIV
       i=1
       do 
          if(i >local_point_factor-1) exit
          N_dim_xyz=N_dim_xyz+(N_dim_cent*3)/2
          N_dim_cent=N_dim_cent*4
          if(i==1) then
             N_dim_xyz_next=N_dim_xyz
             N_dim_cent_next=N_dim_cent
          endif
          i=i+1
       end do
    else
       i=1
       do 
          if(i>local_point_factor-1 .and. & 
               (N_dim_cent>=max_cent .or. do_cavitation) ) exit
          N_dim_xyz=N_dim_xyz+(N_dim_cent*3)/2
          N_dim_cent=N_dim_cent*4
          if(i==1) then
             N_dim_xyz_next=N_dim_xyz
             N_dim_cent_next=N_dim_cent
          endif
          i=i+1
       enddo
    end if
!!!MF <<<<
!     if(local_point_factor-1 < i) local_point_factor=i
!     if(do_cavitation) local_point_factor=1
DPRINT 'maximum allowd centers',N_dim_cent

    allocate(surf_elem%xyz(N_dim_xyz,3),stat=status2)
    ASSERT(status2 == 0)

    !definition of initial set of points
    xyz=>surf_elem%xyz
    xyz=0.0_r8_kind

    xyz(1:4,1)=0.0_r8_kind
    xyz(5,1)=s
    xyz(7,1)=s
    xyz(6,1)=-s
    xyz(8,1)=-s
    xyz(9:10,1)=h
    xyz(11:12,1)=-h
    xyz(13,1)=b 
    xyz(15,1)=b 
    xyz(17,1)=b 
    xyz(19,1)=b 
    xyz(14,1)=-b
    xyz(16,1)=-b
    xyz(18,1)=-b
    xyz(20,1)=-b
    xyz(21,1)=h1
    xyz(22,1)=-h1
    xyz(23,1)=h1
    xyz(24,1)=-h1 
    xyz(25:29,1)=0.0_r8_kind
    xyz(29:30,1)=h2 
    xyz(31:32,1)=-h2
    
    xyz(1,2)=-s
    xyz(3,2)=-s
    xyz(2,2)=s
    xyz(4,2)=s
    xyz(5:6,2)=h
    xyz(7:8,2)=-h
    xyz(9:12,2)=0.0_r8_kind
    xyz(13:14,2)=b
    xyz(17:18,2)=b 
    xyz(15:16,2)=-b 
    xyz(19:20,2)=-b
    xyz(21:24,2)=0.0_r8_kind
    xyz(25,2)=-h2 
    xyz(26:27,2)=h2 
    xyz(28,2)=-h2  
    xyz(29,2)=-h1 
    xyz(30,2)=h1
    xyz(31,2)=-h1 
    xyz(32,2)=h1 
    
    xyz(1:2,3)=h
    xyz(3:4,3)=-h
    xyz(5:8,3)=0.0_r8_kind
    xyz(9,3)=s
    xyz(11,3)=s
    xyz(10,3)=-s
    xyz(12,3)=-s
    xyz(13:16,3)=b 
    xyz(17:20,3)=-b 
    xyz(21:22,3)=h2 
    xyz(23:24,3)=-h2 
    xyz(25:26,3)=h1  
    xyz(27:28,3)=-h1  
    xyz(29:32,3)=0.0_r8_kind
    
    !definition of triangles
    
    allocate(surf_elem%index(N_dim_cent,3),stat=status2)
    ASSERT(status2 == 0)

    indexx=>surf_elem%index
    ind=1
    do i=1,12
       k=i+20
       m=5*(i-1)+1
       indexx(m:m+4,1)=k
       
       n=1
       do j=1,20
          xyz_buf=xyz(k,:)-xyz(j,:)
          r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
          if (r_buf <= 1.06_r8_kind) then
             tang_buf(n)=j
             n=n+1
             if (n > 5) exit
          endif

       enddo

       do l=1,4
          l1=tang_buf(l)
          do n=l+1,5
             n1=tang_buf(n)
             xyz_buf=xyz(n1,:)-xyz(l1,:)
             r_buf=sqrt(dot_product(xyz_buf,xyz_buf))
             if (r_buf <= 2.0_r8_kind*s+0.001_r8_kind) then
                indexx(ind,2)=l1
                indexx(ind,3)=n1
                ind=ind+1
             endif
          enddo
       enddo
    enddo

    ! orientation of points respect to symmetry point group
    if ( (name_point_group=='C1  ') .or. &
         (name_point_group=='C2  ') .or. &
         (name_point_group=='CS  ') .or. &
         (name_point_group=='C2V ') .or. &
         (name_point_group=='C2H ') .or. &
         (name_point_group=='I   ') .or. &
         (name_point_group=='IH  ') .or. &
         (name_point_group=='CI  ')  .or. &
         (name_point_group=='D2  ')  .or. &
         (name_point_group=='D2H ') ) then
       continue !nothing to do
    endif
    if ( (name_point_group=='C3  ') .or. &
         (name_point_group=='C3V ') .or. &
         (name_point_group=='D3  ')  .or. &
         (name_point_group=='D3D ')  .or. &
         (name_point_group=='S6  ')  ) then
       do i=1,N_dim_xyz_st
          x=xyz(i,1)*s/radius+xyz(i,3)*h/radius
          z=-xyz(i,1)*h/radius+xyz(i,3)*s/radius
          xyz(i,1)=x
          xyz(i,3)=z 
       enddo
       if ((name_point_group=='D3  ')  .or. &
            (name_point_group=='D3D ') ) then
          do i=1,N_dim_xyz_st
             x=cos12*xyz(i,1)-sin12*xyz(i,2)
             y=sin12*xyz(i,1)+cos12*xyz(i,2)
             xyz(i,1)=x
             xyz(i,2)=y 
          enddo
       end if
    endif
    if ( (name_point_group=='C5  ') .or. &
         (name_point_group=='C5V ') .or. &
         (name_point_group=='S10 ')  .or. &
         (name_point_group=='D5  ')  .or. &
         (name_point_group=='D5D ')  .or. &
         (name_point_group=='D5H ') ) then
       do i=1,N_dim_xyz_st
          x=xyz(i,1)*cosT+xyz(i,3)*sinT
          z=-xyz(i,1)*sinT+xyz(i,3)*cosT
          xyz(i,1)=x
          xyz(i,3)=z
       enddo
       if ((name_point_group=='D5  ') .or. &
            (name_point_group=='D5D '))  then
          do i=1,N_dim_xyz_st
             x=cos20*xyz(i,1)-sin20*xyz(i,2)
             y=sin20*xyz(i,1)+cos20*xyz(i,2)
             xyz(i,1)=x
             xyz(i,2)=y 
          enddo
       endif
!!!MF <<<<
    endif

!!!MF >>>>
    do j=1,local_point_factor-1
!!!MF allow mor than 1x subdivision
#ifdef _LINUX1
       call more_triangles(N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#else
       call more_triangles(xyz,indexx,N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#endif
       N_dim_xyz_st=N_dim_xyz_next
       N_dim_cent_st=N_dim_cent_next
       N_dim_xyz_next= N_dim_xyz_st+(N_dim_cent_st*3)/2
       N_dim_cent_next= N_dim_cent_st*4
    enddo
!!!MF <<<<

    ! definition of triangle centers
    allocate(surf_elem%xyz_centers(N_dim_cent,3),stat=status2)
    ASSERT(status2 == 0)

    xyz_cent=>surf_elem%xyz_centers
    do i=1,N_dim_cent_st !!!!!!!!N_dim_cent
       i1=indexx(i,1)
       l1=indexx(i,2)
       n1=indexx(i,3)
       xyz_buf=(xyz(i1,:)+xyz(l1,:)+xyz(n1,:))/3.0_r8_kind
       xyz_cent(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
    enddo
    N_points_of_triangles=N_dim_xyz_st !!!!!!N_dim_xyz
    N_centers_on_sphere=N_dim_cent_st  !!!!!!N_dim_cent
    if(do_cavitation .and. output_cavity_data .and. .not.do_gradients) &
         write (output_unit,'(a42,i5)') 'the number of triangles on each sphere is ', N_centers_on_sphere
    
!!$     nullify(xyz,xyz_cent,index)
  end subroutine generate_dodecahedron
  !*************************************************************

  !*************************************************************
!!!MF >>>>
  subroutine generate_doublepyramide(point_factor,do_cavitation,do_gradients)
    !  Purpose: simple extension for all other point groups,
    !           not very well suited, because surface tesserea can have
    !           very acute angles
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),intent(in) :: point_factor
    logical,intent(in)          :: do_cavitation, do_gradients
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind),pointer     :: xyz(:,:),xyz_cent(:,:)
    integer(i4_kind), pointer :: indexx(:,:)
    real(r8_kind)             :: xyz_buf(3)
    integer(i4_kind)          :: N_dim_xyz,N_dim_cent,N_dim_xyz_st,N_dim_cent_st
    integer(i4_kind)          :: N_dim_xyz_next,N_dim_cent_next
    integer(i4_kind)          :: i,j,l1,n1,i1,status3
    real(r8_kind)             :: alph, r_n_rot
    real(r8_kind), parameter  :: pi = 3.14159265355897932368_r8_kind
    real(r8_kind), parameter  :: cos12=0.8660254037844_r8_kind
    real(r8_kind), parameter  :: sin12=0.5_r8_kind
    !------------ Executable code --------------------------------

    radius=1.0_r8_kind 

    N_dim_cent=2*n_rotations
    N_dim_cent_st=2*n_rotations
    N_dim_xyz=n_rotations+2
    N_dim_xyz_st=n_rotations+2

    ! at least one subdivision
    local_point_factor=point_factor+1 !2 is too much !!!!!
    
    i=1
    do 
       if(i>local_point_factor .and. & 
            (N_dim_cent>=max_cent .or. do_cavitation) ) exit
       N_dim_xyz=N_dim_xyz+N_dim_cent*3/2
       N_dim_cent=N_dim_cent*4
       if(i==1) then
          N_dim_xyz_next=N_dim_xyz
          N_dim_cent_next=N_dim_cent
       endif
       i=i+1
    enddo
DPRINT 'maximum allowd centers',N_dim_cent
!!!AS >>>>
!     if(local_point_factor < i) local_point_factor=i-1
!     if(do_cavitation) local_point_factor=0
!!!AS <<<<

    allocate(surf_elem%xyz(N_dim_xyz,3),stat=status3)
    ASSERT(status3 == 0)

    !definition of initial set of points
    xyz=>surf_elem%xyz
    xyz=0.0_r8_kind

    xyz(1,1:2)=0.0_r8_kind
    xyz(1,3)=radius
    xyz(N_dim_xyz_st,1:2)=0.0_r8_kind
    xyz(N_dim_xyz_st,3)=-radius
    xyz(2:N_dim_xyz_st-1,3)=0.0_r8_kind

    if (n_rotations == 3) then
       xyz(2,1)=1.0_r8_kind
       xyz(2,2)=0.0_r8_kind
       xyz(3,1)=-sin12
       xyz(3,2)=cos12
       xyz(4,1)=-sin12
       xyz(4,2)=-cos12
    else
       r_n_rot=real(n_rotations,kind=r8_kind)
       do i=0,n_rotations-1
          alph = (real(i,kind=r8_kind))/r_n_rot * 2.0_r8_kind * pi 
          xyz(i+2,1)=cos(alph)
          xyz(i+2,2)=sin(alph)
       enddo
    endif

    allocate(surf_elem%index(N_dim_cent,3),stat=status3)
    ASSERT(status3 == 0)

    surf_elem%index=0

    indexx=>surf_elem%index
      
    do i=1,n_rotations
       indexx(i,1)=i+1
       indexx(i,2)=i+2
       indexx(i,3)=1
       indexx(n_rotations+i,1)=i+1
       indexx(n_rotations+i,2)=i+2
       indexx(n_rotations+i,3)=N_dim_xyz_st
    enddo
    indexx(n_rotations,2)=2
    indexx(N_dim_cent_st,2)=2
    
    !dividing each triangles in four new triangles 
    do j=1,local_point_factor
!!!MF allow mor than 1x subdivision
#ifdef _LINUX1
       call more_triangles(N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#else
       call more_triangles(xyz,indexx,N_dim_cent_next,N_dim_xyz_st,N_dim_cent_st)
#endif
       N_dim_xyz_st=N_dim_xyz_next
       N_dim_cent_st=N_dim_cent_next
       N_dim_xyz_next= N_dim_xyz_st+N_dim_cent_st*3/2
       N_dim_cent_next= N_dim_cent_st*4
    enddo

    ! definition of triangle centers
    allocate(surf_elem%xyz_centers(N_dim_cent,3),stat=status3)
    ASSERT(status3 == 0)

    xyz_cent=>surf_elem%xyz_centers
    do i=1,N_dim_cent_st
       i1=indexx(i,1)
       l1=indexx(i,2)
       n1=indexx(i,3)
       xyz_buf=(xyz(i1,:)+xyz(l1,:)+xyz(n1,:))/3.0_r8_kind
       xyz_cent(i,:)=(radius/sqrt(dot_product(xyz_buf,xyz_buf)))*xyz_buf
    enddo

    N_points_of_triangles=N_dim_xyz_st
    N_centers_on_sphere=N_dim_cent_st

    if(do_cavitation .and. output_cavity_data .and. .not.do_gradients) &
         write (output_unit,'(a42,i5)') 'the number of triangles on each sphere is ', N_centers_on_sphere
    
!     nullify(xyz,xyz_cent,index)
  end subroutine generate_doublepyramide
!!!MF <<<<
  !*************************************************************

  !*************************************************************
#ifdef _LINUX1
  subroutine more_triangles(n_ind,n_xyz_st,n_ind_st)
#else
  subroutine more_triangles(xyz_t,ind_t,n_ind,n_xyz_st,n_ind_st)
#endif
    !  Purpose: subdivision routine, each triangle is devided in four by
    !           halving the edges
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(r8_kind), pointer          :: xyz_t(:,:)
    integer(i4_kind), pointer       :: ind_t(:,:)
    integer(i4_kind), intent(in)    :: n_ind
    integer(i4_kind), intent(inout) :: n_xyz_st,n_ind_st
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind)                 :: xyz_t_buf(3),xyz_tt(3)
    integer(i4_kind), allocatable :: ind_buf(:,:)
    integer(i4_kind)              :: new_numbers(6)
    integer(i4_kind)              :: status,next_i,neighbour,i1,i2
    integer(i4_kind)              :: i,j,k,l,m,m1,n,n1
    real(r8_kind), parameter      :: small=1.0e-11_r8_kind
    !------------ Executable code --------------------------------

#ifdef _LINUX1
    xyz_t=>surf_elem%xyz
    ind_t=>surf_elem%index
#endif

    allocate(ind_buf(n_ind,3),stat=status)
    ASSERT(status == 0)
    ind_buf=0

    next_i=0
    label_i: do i=1,n_ind_st

       m=0
       label_j: do j=1,3
          next_i=next_i+1
          neighbour=1
          i1=ind_t(i,j)
          ind_buf(next_i,neighbour)=i1

          label_k: do k=1,3
             if (k==j) cycle label_k
             i2=ind_t(i,k)
             xyz_t_buf=(xyz_t(i1,:)+xyz_t(i2,:))/2.0_r8_kind
             xyz_tt=(radius/sqrt(dot_product(xyz_t_buf,xyz_t_buf)))*xyz_t_buf

             label_l: do l=1,n_xyz_st
                
!!$                  if(xyz_tt(1)==xyz_t(l,1) .and. &
!!$                       xyz_tt(2)==xyz_t(l,2) .and. &
!!$                       xyz_tt(3)==xyz_t(l,3)) then
                if(abs(xyz_tt(1)-xyz_t(l,1))<small .and. &
                     abs(xyz_tt(2)-xyz_t(l,2))<small .and. &
                     abs(xyz_tt(3)-xyz_t(l,3))<small) then
                   neighbour=neighbour+1
                   ind_buf(next_i,neighbour)=l
                   m=m+1
                   new_numbers(m)=l
                   cycle label_k
                endif
             enddo label_l
             n_xyz_st=n_xyz_st+1
             xyz_t(n_xyz_st,:)=xyz_tt
             neighbour=neighbour+1
             ind_buf(next_i,neighbour)=n_xyz_st
             m=m+1
             new_numbers(m)=n_xyz_st
          enddo label_k
       enddo label_j
       next_i=next_i+1
       m1=1
       ind_buf(next_i,m1)=new_numbers(1)
       label_n: do n=2,m
          
          label_n1: do n1=1,n-1
             if(new_numbers(n) == new_numbers(n1)) cycle label_n
          enddo label_n1
          m1=m1+1
          ind_buf(next_i,m1)=new_numbers(n)
          if (m1==3) cycle label_i
       enddo label_n
    enddo label_i
    
    n_ind_st=next_i
    ind_t(1:n_ind,1:3)=ind_buf(1:n_ind,1:3)
    
    deallocate(ind_buf,stat=status)
    ASSERT(status == 0)
!     nullify(xyz_t,ind_t)
  end subroutine more_triangles
  !*************************************************************

  !--------------- End of module ----------------------------------
end module polyhedron_module


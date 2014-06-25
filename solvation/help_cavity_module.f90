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
! Public interface of module
!===============================================================
module help_cavity_module
!== Interrupt of public interface of module =========
!  Author: AS
!  Date: 11/05
!
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
!== Interrupt end of public interface of module =====              

!----modules used ------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use group_module ,only: symm_transformation_int,ylm_trafos,sub_group, &
       group_coset,group_num_el,group_coset_decomp
  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------
  public calc_symm_pairs,calc_intersect_spheres,if_tess_inside_cav,center_inside_cav

!== Interrupt of public interface of module =========
  real(r8_kind) , parameter :: pi = 3.14159265355897932368_r8_kind

  external error_handler

contains

  !********************************************************************
  subroutine calc_symm_pairs(in1, in2, drr, r, NN, pairs, NP, rns, r_sym)
    !
    ! FIXME: literal constant (500) in interface definition?
    !
    implicit none
    integer(i4_kind), intent(in) :: in1, in2, NN
    integer(i4_kind), intent(out) :: NP
    integer(i4_kind) :: pairs(500, 2) ! FIXME: intent(out)?
    real(r8_kind), intent(in) :: r(NN, 3) ! (N_spheres, 3)
    real(r8_kind), intent(in) :: drr
    real(r8_kind), intent(in), optional :: rns(3)
    real(r8_kind), optional :: r_sym(500, 4) ! FIXME: intent(out)?
    ! *** end of interface ***

    !------------ Declaration of local variables ------------------
    integer(i4_kind) :: pairs_tmp(500,2)
    real(r8_kind)    :: r_sym_tmp(500,3)
    real(r8_kind)    :: position1_1(3),position2_1(3),position3_1(3)
    real(r8_kind)    :: position1_2(3),position2_2(3),position3_2(3)
    real(r8_kind)    :: position1_3(3),position2_3(3),position3_3(3)
    type(sub_group) :: local_groups_1,local_groups_2,local_groups_3
    type(group_coset),target :: cosets_1,cosets_2,cosets_3
    type(symm_transformation_int) :: point_trafos
    real(r8_kind)    :: rn(3),dr1,dr2,xyz(3),xyz1(3),xyz2(3)
    real(r8_kind)    :: d1,d2,dd

    integer(i4_kind) :: n_1,n_2,n_3,n_equal_1,n_equal_2,n_equal_3
    integer(i4_kind) :: i,j,k,l,m,n,status
    real(r8_kind),parameter :: small = 1.0e-10_r8_kind
    !------------ Executable code ---------------------------------

    pairs_tmp=0
    r_sym_tmp=0.0_r8_kind
    ! reorder coordinates of points as
    ! (x,y,z) --> (x,z,y) in order to comply with the
    ! convention for angular momentum l=1
!!$print*,r(in1,:),'1'
!!$print*,r(in2,:),'2'

    position1_1(1) = r(in1,1); position1_1(2) = r(in1,3); position1_1(3) = r(in1,2)
    position1_2(1) = r(in2,1); position1_2(2) = r(in2,3); position1_2(3) = r(in2,2)
    if(.not.present(rns)) then 
       rn=(r(in1,:)+r(in2,:))/2.0_r8_kind
    else
       rn=rns
    end if
!!$print*,rn,'3'
    position1_3(1) = rn(1); position1_3(2) = rn(3); position1_3(3) = rn(2)
    dr1=sqrt(dot_product(rn-r(in1,:),rn-r(in1,:)))
    dr2=sqrt(dot_product(rn-r(in2,:),rn-r(in2,:)))

    ! determine local symmetry groups
    !
    ! now apply all symmetry operations to the position of the points
    n_1=0; n_2=0; n_3=0
    do j=1,group_num_el
       position2_1=MATMUL(ylm_trafos(1)%matrix(:,:,j),position1_1)
       if (dot_product(position2_1-position1_1,position2_1-position1_1) < small) then
          n_1=n_1+1
       endif
       position2_2=MATMUL(ylm_trafos(1)%matrix(:,:,j),position1_2)
       if (dot_product(position2_2-position1_2,position2_2-position1_2) < small) then
          n_2=n_2+1
       endif
       position2_3=MATMUL(ylm_trafos(1)%matrix(:,:,j),position1_3)
       if (dot_product(position2_3-position1_3,position2_3-position1_3) < small) then
          n_3=n_3+1
       endif
    enddo

    ! allocate group elements
    local_groups_1%num_el=n_1; local_groups_2%num_el=n_2; local_groups_3%num_el=n_3
    allocate(local_groups_1%elements(local_groups_1%num_el), &
         local_groups_2%elements(local_groups_2%num_el), &
         local_groups_3%elements(local_groups_3%num_el),stat=status)
    if ( status /= 0) call error_handler( &
         "help_cavity_module:calc_symm_pairs: failed allocation of local_groups")

    ! fill up group elements
    n_1=0
    do j=1,group_num_el
       position2_1=MATMUL(ylm_trafos(1)%matrix(:,:,j),position1_1)
       if (dot_product(position2_1-position1_1,position2_1-position1_1) < small) then
          n_1=n_1+1
          local_groups_1%elements(n_1)=j
       end if
    enddo
    n_2=0
    do j=1,group_num_el
       position2_2=MATMUL(ylm_trafos(1)%matrix(:,:,j),position1_2)
       if (dot_product(position2_2-position1_2,position2_2-position1_2) < small) then
          n_2=n_2+1
          local_groups_2%elements(n_2)=j
       end if
    enddo
    n_3=0
    do j=1,group_num_el
       position2_3=MATMUL(ylm_trafos(1)%matrix(:,:,j),position1_3)
       if (dot_product(position2_3-position1_3,position2_3-position1_3) < small) then
          n_3=n_3+1
          local_groups_3%elements(n_3)=j
       end if
    enddo

    ! now determine symmetry equivalent atoms
    call group_coset_decomp(n_equal_1,local_groups_1,cosets_1,point_trafos%matrix)
    call group_coset_decomp(n_equal_2,local_groups_2,cosets_2,point_trafos%matrix)
    call group_coset_decomp(n_equal_3,local_groups_3,cosets_3,point_trafos%matrix)

    !calculate symmetry equivalent pairs
    j_n_equal_1: do j=1,n_equal_1
       position2_1=MATMUL(ylm_trafos(1)%matrix(:,:,cosets_1%elements(1,j)),position1_1)
       position3_1(1)=position2_1(1)
       position3_1(2)=position2_1(3)
       position3_1(3)=position2_1(2)
       do i=1,NN
          if(sqrt(dot_product(position3_1-r(i,:),position3_1-r(i,:))) <= small) exit
       end do
!!$print*,position3_1,'11',i
       pairs_tmp(j,1)=i
    end do j_n_equal_1

    j_n_equal_2: do j=1,n_equal_2
       position2_2=MATMUL(ylm_trafos(1)%matrix(:,:,cosets_2%elements(1,j)),position1_2)
       position3_2(1)=position2_2(1)
       position3_2(2)=position2_2(3)
       position3_2(3)=position2_2(2)
       do i=1,NN
          if(sqrt(dot_product(position3_2-r(i,:),position3_2-r(i,:))) <= small) exit
       end do
!!$print*,position3_2,'22',i
       pairs_tmp(j,2)=i
    end do j_n_equal_2

    j_n_equal_3: do j=1,n_equal_3
       position2_3=MATMUL(ylm_trafos(1)%matrix(:,:,cosets_3%elements(1,j)),position1_3)
       position3_3(1)=position2_3(1)
       position3_3(2)=position2_3(3)
       position3_3(3)=position2_3(2)
!!$print*,position3_3,'33',j
       r_sym_tmp(j,:)=position3_3
    end do j_n_equal_3

    deallocate(local_groups_1%elements, &
         local_groups_2%elements, &
         local_groups_3%elements,stat=status)
    if ( status /= 0) call error_handler( &
         "help_cavity_module:calc_symm_pairs: failed deallocation of local_groups")

    n=0
    do i=1,n_equal_3
       xyz=r_sym_tmp(i,:)
       do j=1,n_equal_1
          l=pairs_tmp(j,1)
          xyz1=r(l,:)
          d1=sqrt(dot_product(xyz-xyz1,xyz-xyz1))
          if(abs(dr1-d1) <= small) then
             do k=1,n_equal_2
                m=pairs_tmp(k,2)
                if(m==l) cycle
                xyz2=r(m,:)
                d2=sqrt(dot_product(xyz-xyz2,xyz-xyz2))
                dd=sqrt(dot_product(xyz1-xyz2,xyz1-xyz2))
                if((abs(dr2-d2) <= small) .and. (abs(drr-dd) <= small)) then
                   n=n+1
                   pairs(n,1)=l; pairs(n,2)=m
                end if
             end do
          end if
       end do
    end do

    NP=n
    do i=1,NP
       if(pairs(i,2) < pairs(i,1)) then
          k=pairs(i,1); pairs(i,1)=pairs(i,2); pairs(i,2)=k
       end if
    end do

    pairs_tmp=pairs; pairs=0

    pairs(1,:)=pairs_tmp(1,:)
    k=1
    c2:do i=2,NP
       do j=1,i-1
          if((pairs_tmp(i,1)==pairs(j,1) .and. pairs_tmp(i,2)==pairs(j,2))) cycle c2
       end do
       k=k+1
       pairs(k,:)=pairs_tmp(i,:)
    end do c2
    
    NP=k

    do i=1,NP
!!$print*,pairs(i,:)
       if(present(rns)) then
          r_sym(i,1:3)=r_sym_tmp(i,:)
          r_sym(i,4)=1.0_r8_kind
       end if
    end do

    if(present(rns)) then
       l1:do i=2,NP
          do j=1,i-1
             if(sqrt(dot_product(r_sym(i,1:3)-r_sym(j,1:3), &
                  r_sym(i,1:3)-r_sym(j,1:3))) <= small) then
                r_sym(i,4)=0.0_r8_kind
                cycle l1
             end if
          end do
       end do l1
    end if

  end subroutine calc_symm_pairs
  !********************************************************************

  !********************************************************************
  subroutine calc_intersect_spheres(i_sphere,xyz_spheres,r_spheres,NN,n_spheres,j_spheres)

    integer(i4_kind) :: i_sphere,NN,n_spheres
    real(r8_kind) :: xyz_spheres(NN,3),r_spheres(NN)
    integer(i4_kind) :: j_spheres(NN)
    !------------ Declaration of local variables ------------------
    integer(i4_kind) :: i,j,k
    real(r8_kind) :: dist
    !------------ Executable code ---------------------------------

    k=0
    do i=1,NN
       j=i_sphere
       if(i==j) cycle
       dist=sqrt(dot_product(xyz_spheres(j,:)-xyz_spheres(i,:),xyz_spheres(j,:)-xyz_spheres(i,:)))
       if(dist <= r_spheres(j)+r_spheres(i)) then
          k=k+1
          j_spheres(k)=i
       end if
    end do

    n_spheres=k

  end subroutine calc_intersect_spheres
  !********************************************************************

  !********************************************************************
  function center_inside_cav(s,xyz_spheres,r_spheres,NN,n_spheres,j_spheres)

    logical :: center_inside_cav

    real(r8_kind) :: s(3)
    integer(i4_kind) :: NN,n_spheres
    real(r8_kind) :: xyz_spheres(NN,3),r_spheres(NN)
    integer(i4_kind) :: j_spheres(NN)
    !------------ Declaration of local variables ------------------
    real(r8_kind) :: dist
    logical :: outside
    real(r8_kind),parameter :: small = 1.0e-8_r8_kind
    integer(i4_kind) :: i,j
    !------------ Executable code ---------------------------------

    if(n_spheres==0) then
       center_inside_cav=.false.
       return
    end if

    outside=.false.
    do i=1,n_spheres
       j=j_spheres(i)
       dist=sqrt(dot_product(s(:)-xyz_spheres(j,:),s(:)-xyz_spheres(j,:)))
       if(dist <= r_spheres(j)+small) then
          outside=.false.
          exit
       end if
       outside=.true.
    end do
    if(outside) then
       center_inside_cav=.false.
    else
       center_inside_cav=.true.
    end if

  end function center_inside_cav
  !********************************************************************

  !********************************************************************
  function if_tess_inside_cav(v1, v2, v3, i_sphere, N_t, xyz_spheres, r_spheres, NN, n_spheres, j_spheres)
    implicit none
    logical :: if_tess_inside_cav
    real(r8_kind) :: v1(:), v2(:), v3(:) ! (3)
    integer(i4_kind) :: i_sphere,N_t,NN,n_spheres
    real(r8_kind) :: xyz_spheres(NN,3),r_spheres(NN)
    integer(i4_kind) :: j_spheres(NN)
    ! *** end of interface ***

    !------------ Declaration of local variables ------------------
    real(r8_kind) :: dist,area
    integer(i4_kind),parameter ::n_tessel=4
    logical :: outside
    real(r8_kind) :: sub_tes(256,3,3),sub_tes_buf(256,3,3),xyz(3,3),cent(3)
    real(r8_kind),parameter :: small = 1.0e-8_r8_kind
    integer(i4_kind) :: i,j,k
    !------------ Executable code ---------------------------------

    !Check if the initial tessarea has all its 3 vertex within cavity.
    !if not cancel checking to begin tesselation procedure
    if(n_spheres==0) then
       if_tess_inside_cav=.false.
       return
    end if

    outside=.false.
    do i=1,n_spheres
       j=j_spheres(i)
       dist=sqrt(dot_product(v1(:)-xyz_spheres(j,:),v1(:)-xyz_spheres(j,:)))
       if(dist <= r_spheres(j)+small) then
          outside=.false.
          cycle
       end if
       outside=.true.
    end do
    if(outside) then
       if_tess_inside_cav=.false.
       return
    end if

    outside=.false.
    do i=1,n_spheres
       j=j_spheres(i)
       dist=sqrt(dot_product(v2(:)-xyz_spheres(j,:),v2(:)-xyz_spheres(j,:)))
       if(dist <= r_spheres(j)+small) then
          outside=.false.
          cycle
       end if
       outside=.true.
    end do
    if(outside) then
       if_tess_inside_cav=.false.
       return
    end if

    outside=.false.
    do i=1,n_spheres
       j=j_spheres(i)
       dist=sqrt(dot_product(v3(:)-xyz_spheres(j,:),v3(:)-xyz_spheres(j,:)))
       if(dist <= r_spheres(j)+small) then
          outside=.false.
          cycle
       end if
       outside=.true.
    end do
    if(outside) then
       if_tess_inside_cav=.false.
       return
    end if

    area=(4.0_r8_kind*pi*r_spheres(i_sphere)**2)/N_t

    !calculating grid on tessarea surface
    
    !calculating subtessareas
    sub_tes(1,1,:)=v1
    sub_tes(1,2,:)=v2
    sub_tes(1,3,:)=v3
    do i=1,n_tessel
       k=0
       do j=1,4**(i-1)
          xyz(1,:)=(sub_tes(j,1,:)+sub_tes(j,2,:))/2.0_r8_kind
          xyz(1,:)=(r_spheres(i_sphere)/sqrt(dot_product(xyz(1,:),xyz(1,:))))*xyz(1,:)
          xyz(2,:)=(sub_tes(j,1,:)+sub_tes(j,3,:))/2.0_r8_kind
          xyz(2,:)=(r_spheres(i_sphere)/sqrt(dot_product(xyz(2,:),xyz(2,:))))*xyz(2,:)
          xyz(3,:)=(sub_tes(j,2,:)+sub_tes(j,3,:))/2.0_r8_kind
          xyz(3,:)=(r_spheres(i_sphere)/sqrt(dot_product(xyz(3,:),xyz(3,:))))*xyz(3,:)

          k=k+1
          sub_tes_buf(k,1,:)=sub_tes(j,1,:); sub_tes_buf(k,2,:)=xyz(1,:); sub_tes_buf(k,3,:)=xyz(2,:)
          k=k+1
          sub_tes_buf(k,1,:)=sub_tes(j,2,:); sub_tes_buf(k,2,:)=xyz(1,:); sub_tes_buf(k,3,:)=xyz(3,:)
          k=k+1
          sub_tes_buf(k,1,:)=sub_tes(j,3,:); sub_tes_buf(k,2,:)=xyz(2,:); sub_tes_buf(k,3,:)=xyz(3,:)
          k=k+1
          sub_tes_buf(k,1,:)=xyz(1,:); sub_tes_buf(k,2,:)=xyz(2,:); sub_tes_buf(k,3,:)=xyz(3,:)
       end do
       sub_tes(1:k,:,:)=sub_tes_buf(1:k,:,:)
    end do

    !calculating grid on tessarea surface
    outside=.false.
    do i=1,4**n_tessel
       cent(:)=(sub_tes(i,1,:)+sub_tes(i,2,:)+sub_tes(i,3,:))/3.0_r8_kind
       cent(:)=(r_spheres(i_sphere)/sqrt(dot_product(cent(:),cent(:))))*cent(:)
       do j=1,n_spheres
          k=j_spheres(j)
          dist=sqrt(dot_product(cent(:)-xyz_spheres(k,:),cent(:)-xyz_spheres(k,:)))
          if(dist <= r_spheres(j)+small) then
             outside=.false.
             exit
          end if
          outside=.true.
       end do
       if(outside) exit
    end do
    
    if(outside) then
       if_tess_inside_cav=.false.
    else
       if_tess_inside_cav=.true.
    end if

  end function if_tess_inside_cav
  !********************************************************************

  !********************************************************************

end module help_cavity_module


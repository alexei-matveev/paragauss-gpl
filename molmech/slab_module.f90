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
module slab_module
  !------------ Modules used -----------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use tasks_main_options_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module ====================
  !------------ Declaration of public constants and variables -----
  type, public :: slab_cell
     real(kind=r8_kind) :: a, b
     real(kind=r8_kind) :: alpha
  end type slab_cell

  type, public :: slab_vectors
     real(kind=r8_kind) :: v1(2), v2(2)
  end type slab_vectors

  type, public :: slab_image_coor
     integer(i4_kind) :: type
     real(r8_kind), pointer :: r(:,:)
  end type slab_image_coor

  type(slab_cell), public :: cel_s
  type(slab_vectors), public :: vect_s  !vectors of real slab unit cell
  type(slab_vectors), public :: kvect_s !vectors of reciprocal slab unit cell
  real(r8_kind), public :: area
  type(slab_image_coor), public, allocatable :: im_coor_s(:)
  integer(kind=i4_kind), public :: n_images_s

  real(r8_kind), public :: strain_s(3) ! 1-xx,2-yy, 3-xy

  logical, public :: slab_calc
  !------------ public functions and subroutines ---------------------
  public read_slab_cell, write_slab_cell_vect_to_output, read_slab_vectors, &
       cart2frac_slab, frac2cart_slab, valid_cutoff_slab, image_slab, &
       renew_slab, store_vect_s, use_vect_s
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of private constants and variables ----
  type(slab_vectors) :: vect_s_store  !store buffer for vectors of real unit cell

  real(kind=r8_kind) :: a, b, alpha 

  namelist /cell_slab/ a, b, alpha

  real(kind=r8_kind) :: vector1(2), vector2(2)
  
  namelist /vectors_slab/ vector1, vector2

  !------------ Subroutines ------------------------------------------
contains
  !******************************************************************
  function read_slab_cell()

    logical :: read_slab_cell
    integer(i4_kind) :: i

    call go_to_first_input_line
    read_slab_cell=find_namelist("&CELL_SLAB",i)
    if(.not.read_slab_cell) return
    read(input_device,nml=cell_slab, end=100, err=200)

    if(a <= zero .or. b <= zero .or. alpha <= zero) goto 200
    cel_s%a=a; cel_s%b=b
    cel_s%alpha=alpha
    read_slab_cell=.true.

    call cell_s2vect_s(cel_s%a,cel_s%b,cel_s%alpha,vect_s%v1,vect_s%v2)

    area=slab_area()
    call calc_k_vectors_s(vect_s%v1,vect_s%v2,area,kvect_s%v1,kvect_s%v2)
    return

100 read_slab_cell=.false.
    return
200 call input_nm_error(0,"CELL_SLAB") 

  end function read_slab_cell
  !******************************************************************

  !******************************************************************
  function read_slab_vectors()

    logical :: read_slab_vectors
    real(kind=r8_kind) :: r1
    integer(i4_kind) :: i

    call go_to_first_input_line
    read_slab_vectors=find_namelist("&VECTORS_SLAB",i)
    if(.not.read_slab_vectors) return
    read(input_device,nml=vectors_slab, end=100, err=200)

    r1=sqrt(dot_product(vector1,vector1))
    if(r1 == zero) goto 200
    vect_s%v1=vector1; vect_s%v2=vector2

    call vect_s2cell_s(vect_s%v1,vect_s%v2,cel_s%a,cel_s%b,cel_s%alpha)
    read_slab_vectors=.true.

    area=slab_area()
    call calc_k_vectors_s(vect_s%v1,vect_s%v2,area,kvect_s%v1,kvect_s%v2)
    return
    
100 read_slab_vectors=.false.
    return
200 call input_nm_error(0,"VECTORS_SLAB") 

  end function read_slab_vectors
  !******************************************************************

  !******************************************************************
  subroutine write_slab_cell_vect_to_output()

    write(output_device,'(a29,/)') " Parameters of the SLAB unit cell:"
    write(output_device,'(2(a5,f13.8),a11,/)')  &
         "a = ",cel_s%a,"b= ",cel_s%b,"(angstrom)"
    write(output_device,'(a8,f12.8,a9,/)')  &
         "alpha= ",cel_s%alpha,"(degree)"

    
    write(output_device,'(a43,/)') &
         " Cartesian coordinates of the slab vectors:"
    write(output_device,'(2(1x,f13.8))') vect_s%v1
    write(output_device,'(2(1x,f13.8),/)') vect_s%v2

    write(output_device,'(a17,f20.10,a11,/)') ' Unit cell area: ',area,' angstrom^2'

    if(minimal_image) write(output_device,'(a65,/)') & 
         ' Periodic boundary conditions: the minimal image technics is used'

  end subroutine write_slab_cell_vect_to_output
  !******************************************************************

  !******************************************************************
  subroutine cell_s2vect_s(ac,bc,alphac,vc1,vc2)

    real(kind=r8_kind), intent(in) :: ac,bc,alphac
    real(kind=r8_kind), intent(out) :: vc1(2),vc2(2)

    real(kind=r8_kind) :: alpha_r

    alpha_r=pi*alphac/pi_degree

    vc1(1)=ac; vc1(2)=zero
    vc2(1)=bc*cos(alpha_r); vc2(2)=bc*sin(alpha_r)

  end subroutine cell_s2vect_s
  !******************************************************************

  !******************************************************************
  subroutine vect_s2cell_s(vc1,vc2,ac,bc,alphac)

    real(kind=r8_kind), intent(in) :: vc1(2),vc2(2)
    real(kind=r8_kind), intent(out) :: ac,bc,alphac

    real(kind=r8_kind) :: alpha_r
    integer(kind=i4_kind) :: i

    ac=sqrt(dot_product(vc1,vc1))
    bc=sqrt(dot_product(vc2,vc2))

    alpha_r=zero
    do i=1,2
       alpha_r=alpha_r+vc2(i)*vc1(i)
    end do

    alpha_r=acos(alpha_r/(ac*bc))

    alphac=pi_degree*alpha_r/pi
    
  end subroutine vect_s2cell_s
  !******************************************************************

  !******************************************************************
  function slab_area()

    real(r8_kind) :: slab_area
 
    slab_area=vect_s%v1(1)*vect_s%v2(2)-vect_s%v2(1)*vect_s%v1(2)

  end function slab_area
  !******************************************************************

  !******************************************************************
  subroutine calc_k_vectors_s(vc1,vc2,s,kv1,kv2)

    real(kind=r8_kind), intent(in) :: vc1(2),vc2(2),s
    real(kind=r8_kind), intent(out) :: kv1(2),kv2(2)
    
    kv1(1)=two*pi*vc2(2)/s; kv1(2)=-two*pi*vc2(1)/s
    kv2(1)=-two*pi*vc1(2)/s; kv2(2)=two*pi*vc1(1)/s

  end subroutine calc_k_vectors_s
  !******************************************************************

  !******************************************************************
  function valid_cutoff_slab(rcut)

    real(r8_kind) :: valid_cutoff_slab
    real(r8_kind),intent(in) :: rcut
    real(r8_kind) :: rmax
    real(r8_kind) :: amax,bmax,a_sin

    a_sin=sin(pi*cel_s%alpha/pi_degree)
    
    amax=cel_s%a*a_sin
    bmax=cel_s%b*a_sin

    rmax=half*min(amax,bmax)

    if(rcut <= rmax) then
       valid_cutoff_slab=zero
    else
       valid_cutoff_slab=rmax
    end if

  end function valid_cutoff_slab
  !******************************************************************

  !******************************************************************
  subroutine image_slab(r)
    
    real(r8_kind), intent(inout) :: r(3)

    real(r8_kind) :: fr(3),bmat(2,2),det,rr
    integer(i4_kind) :: i

    bmat(1,1)=vect_s%v2(2)
    bmat(2,1)=-vect_s%v1(2)
    bmat(1,2)=-vect_s%v2(1)
    bmat(2,2)=vect_s%v1(1)

    det=vect_s%v1(1)*vect_s%v2(2)-vect_s%v2(1)*vect_s%v1(2)
    rr=zero
    if(abs(det) > zero) rr=one/det

    bmat=bmat*rr

    fr(1:2)=matmul(bmat,r(1:2))
    fr(3)=r(3)
    do i=1,2
       if(abs(fr(i)) <= half-small1 .or. abs(fr(i)) >= half+small1) then 
          fr(i)=fr(i)-anint(fr(i))
       end if
    end do
    
    bmat(:,1)=vect_s%v1
    bmat(:,2)=vect_s%v2

    r(1:2)=matmul(bmat,fr(1:2))
   
  end subroutine image_slab
  !******************************************************************

  !******************************************************************
  subroutine renew_slab(new_strain)

    real(r8_kind) :: new_strain(3)
    real(r8_kind) :: strain_mat(2,2)
    real(r8_kind) :: v1(2),v2(2),v1p(2),v2p(2)
    integer(i4_kind) :: i

    strain_mat(1,1)=one+new_strain(1)
    strain_mat(2,2)=one+new_strain(2)
    strain_mat(1,2)=half*new_strain(3)
    strain_mat(2,1)=strain_mat(1,2)

    v1p=vect_s%v1; v2p=vect_s%v2

    v1=zero; v2=zero
    do i=1,3
       v1(i)=v1(i)+strain_mat(1,i)*v1p(1)+strain_mat(2,i)*v1p(2)
       v2(i)=v2(i)+strain_mat(1,i)*v2p(1)+strain_mat(2,i)*v2p(2)
    end do

    vect_s%v1=v1; vect_s%v2=v2

    call vect_s2cell_s(vect_s%v1,vect_s%v2, &
         cel_s%a,cel_s%b,cel_s%alpha)
    area=slab_area()
    call calc_k_vectors_s(vect_s%v1,vect_s%v2,area,kvect_s%v1,kvect_s%v2)

  end subroutine renew_slab
  !******************************************************************

  !******************************************************************
  subroutine store_vect_s()

    vect_s_store=vect_s

  end subroutine store_vect_s
  !******************************************************************

  !******************************************************************
  subroutine use_vect_s()

    vect_s=vect_s_store

  end subroutine use_vect_s
  !******************************************************************

  !******************************************************************
  subroutine cart2frac_slab(cart,frac)

    real(r8_kind) :: cart(3),frac(3)

    real(kind=r8_kind) :: bmat(2,2),det,rr
    integer(kind=i4_kind) :: i

    bmat(1,1)=vect_s%v2(2)
    bmat(2,1)=-vect_s%v1(2)
    bmat(1,2)=-vect_s%v2(1)
    bmat(2,2)=vect_s%v1(1)

    det=vect_s%v1(1)*vect_s%v2(2)-vect_s%v2(1)*vect_s%v1(2)
    rr=zero
    if(abs(det) > zero) rr=one/det

    bmat=bmat*rr

    frac(1:2)=matmul(bmat,cart(1:2))
    frac(3)=cart(3)
    do i=1,2
       if(abs(frac(i)) <= half-small1 .or. abs(frac(i)) >= half+small1) then
          frac(i)=frac(i)-anint(frac(i))
       end if
       if(abs(frac(i)) <= small1) frac(i)=zero
       if(frac(i) < zero) then 
          frac(i)=one+frac(i)
       elseif(frac(i) >= one) then
          frac(i)=frac(i)-one
       end if
    end do

  end subroutine cart2frac_slab
  !******************************************************************

  !******************************************************************
  subroutine frac2cart_slab(frac,cart)

    real(r8_kind) :: cart(3),frac(3)

    real(kind=r8_kind) :: bmat(2,2)

    bmat(:,1)=vect_s%v1
    bmat(:,2)=vect_s%v2

    cart(1:2)=matmul(bmat,frac(1:2))
    cart(3)=frac(3)

  end subroutine frac2cart_slab
  !******************************************************************
end module slab_module










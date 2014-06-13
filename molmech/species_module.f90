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
module species_module
# include <def.h>
  !------------ Modules used --------------------------------------
  use type_module
  use common_data_module
  use inp_out_module
  use element_data_module
  use tasks_main_options_module
  use slab_module
  use qmmm_interface_module
  use molmech_msgtag_module
  use comm_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  type, public :: atom_coord
     integer(kind=i4_kind) :: type
     integer(kind=i4_kind) :: initial_number
     real(kind=r8_kind) :: r(3)
     integer(kind=i4_kind) :: neigh_sp
  end type atom_coord
  
  type, public :: at_species
     character(len=len_name) :: name
     integer(kind=i4_kind) :: main_number
     integer(kind=i4_kind) :: c_s
     real(kind=r8_kind) :: charge
     real(kind=r8_kind) :: r_coval
     real(kind=r8_kind) :: r_vdw
     real(kind=r8_kind) :: redfac
     real(kind=r8_kind) :: r0_dr
     real(kind=r8_kind) :: k_dr
     real(kind=r8_kind) :: scalefac
     real(kind=r8_kind) :: scalefac_smooth
     logical :: no_sphere
  end type at_species

  type, public :: lat_cell
     real(kind=r8_kind) :: a, b, c
     real(kind=r8_kind) :: alpha, beta, gamma
  end type lat_cell

  type, public :: lat_vectors
     real(kind=r8_kind) :: v1(3), v2(3), v3(3)
  end type lat_vectors

  type, public :: image_coor
     integer(i4_kind) :: type
     real(r8_kind), pointer :: r(:,:)
  end type image_coor

  type(atom_coord), allocatable, public :: atoms_cart(:)
  type(atom_coord), allocatable, public :: atoms_frac(:)
  type(at_species), allocatable, public :: atoms(:)

  type(lat_cell), public :: cel
  type(lat_vectors), public :: vect  !vectors of real unit cell
  type(lat_vectors), public :: kvect !vectors of reciprocal unit cell
  real(r8_kind), public :: volume
  type(image_coor), public, allocatable :: im_coor(:)
  type(image_coor), public, allocatable :: im_coor_solv(:)
  integer(kind=i4_kind), public :: n_images,n_images_solv
  real(r8_kind), public :: strain(6) !1-xx,2-yy,3-zz,6-xy,5-xz,4-yz

  integer(i4_kind), allocatable, public :: cs_pair(:,:) !to prepare gx.cv file (epe interface)
  ! definition in n_body_lists_module

  integer(kind=i4_kind), public :: n_species, n_species_types
  logical, public :: lattice_calc
  integer(kind=i4_kind), public :: resort_axis(1)
  real(r8_kind), public :: gx_work
  integer, public :: n_fixed
  !------------ public functions and subroutines ------------------
  public read_coord, read_coord_gx, read_species, read_pdc,write_species_to_output, &
       name2type, get_c_s, read_cell, read_vectors, write_cell_vect_to_output, &
       cart2frac, frac2cart, species_resorting, write_final_geometry, save_xyz, save_vm, &
       epe_interface, read_qmmm, shutdown_species_data, valid_cutoff, &
       image, calc_images_species, calc_images_slab_species, send_receive_species, &
       renew_lattice, store_vect, use_vect, check_total_charge_and_dipmom, &
       shutdown_species_on_slaves, calc_images_slab_species_solvation, &
       calc_images_species_solvation
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  type(lat_vectors) :: vect_store  !store buffer for vectors of real unit cell

  real(kind=r8_kind) :: r_coval
  real(kind=r8_kind) :: r_vdw
  real(kind=r8_kind) :: r0_dr
  real(kind=r8_kind) :: k_dr
  real(kind=r8_kind) :: charge
  character(len=len_name) :: name
  character(len=2) :: main_name
  real(kind=r8_kind) :: vdw_red_fac
  real(kind=r8_kind) :: solv_scale_fac,solv_scale_fac_smooth
  logical :: no_sphere

  real(kind=r8_kind) :: df_charge=0.0_r8_kind
  real(kind=r8_kind) :: df_r_coval=0.0_r8_kind
  real(kind=r8_kind) :: df_r_vdw=0.0_r8_kind
  real(kind=r8_kind) :: df_r0_dr=zero
  real(kind=r8_kind) :: df_k_dr=zero
  character(len=len_name) :: df_name="    "
  character(len=2) :: df_main_name="XX"
  real(kind=r8_kind) :: df_vdw_red_fac=one
  real(kind=r8_kind) :: df_solv_scale_fac=1.2_r8_kind
  real(kind=r8_kind) :: df_solv_scale_fac_smooth=1.125_r8_kind
  logical :: df_no_sphere=.false.

  namelist /species/ charge,name,r_coval,r_vdw,main_name,vdw_red_fac,r0_dr,k_dr, &
       solv_scale_fac,solv_scale_fac_smooth,no_sphere

  real(kind=r8_kind) :: a, b, c, alpha, beta, gamma

  namelist /cell/ a, b, c, alpha, beta, gamma

  real(kind=r8_kind) :: vector1(3), vector2(3), vector3(3)
  
  namelist /vectors/ vector1, vector2, vector3

  !------------ Subroutines ---------------------------------------
contains
  !******************************************************************
  function read_cell()

    logical :: read_cell
    integer(i4_kind) :: i

    call go_to_first_input_line
    read_cell=find_namelist("&CELL",i)
    if(.not.read_cell) return
    read(input_device,nml=cell, end=100, err=200)

    if(a <= zero .or. b <= zero .or. c <= zero .or. &
         alpha <= zero .or. beta <= zero .or. gamma <= zero) goto 200
    cel%a=a; cel%b=b; cel%c=c
    cel%alpha=alpha; cel%beta=beta; cel%gamma=gamma
    read_cell=.true.

    call cell2vect(cel%a,cel%b,cel%c,cel%alpha,cel%beta,cel%gamma, &
         vect%v1,vect%v2,vect%v3)

    volume=cell_volume()
    call calc_k_vectors(vect%v1,vect%v2,vect%v3,volume,kvect%v1,kvect%v2,kvect%v3)
    return

100 read_cell=.false.
    return
200 call input_nm_error(0,"CELL") 

  end function read_cell
  !******************************************************************

  !******************************************************************
  function read_vectors()

    logical :: read_vectors
    real(kind=r8_kind) :: r1,r2,r3
    integer(i4_kind) :: i

    call go_to_first_input_line
    read_vectors=find_namelist("&VECTORS",i)
    if(.not.read_vectors) return
    read(input_device,nml=vectors, end=100, err=200)

    r1=sqrt(dot_product(vector1,vector1))
    r2=sqrt(dot_product(vector2,vector2))
    r3=sqrt(dot_product(vector3,vector3))
    if(r1 == zero .or. r2 == zero .or. r3 == zero) goto 200
    vect%v1=vector1; vect%v2=vector2; vect%v3=vector3

    call vect2cell(vect%v1,vect%v2,vect%v3, &
         cel%a,cel%b,cel%c,cel%alpha,cel%beta,cel%gamma)
    read_vectors=.true.

    volume=cell_volume()
    call calc_k_vectors(vect%v1,vect%v2,vect%v3,volume,kvect%v1,kvect%v2,kvect%v3)
    return
    
100 read_vectors=.false.
    return
200 call input_nm_error(0,"VECTORS") 

  end function read_vectors
  !******************************************************************

  !******************************************************************
  subroutine write_cell_vect_to_output()

    write(output_device,'(a29,/)') " Parameters of the unit cell:"
    write(output_device,'(3(a5,f13.8),a11,/)')  &
         "a= ",cel%a,"b= ",cel%b,"c= ",cel%c,"(angstrom)"
    write(output_device,'(3(a8,f12.8),a9,/)')  &
         "alpha= ",cel%alpha,"beta = ",cel%beta,"gamma= ",cel%gamma,"(degree)"

    
    write(output_device,'(a46,/)') &
         " Cartesian coordinates of the lattice vectors:"
    write(output_device,'(3(1x,f13.8))') vect%v1
    write(output_device,'(3(1x,f13.8))') vect%v2
    write(output_device,'(3(1x,f13.8),/)') vect%v3

    write(output_device,'(a19,f20.10,a11,/)') ' Unit cell volume: ',volume,' angstrom^3'

    if(minimal_image) write(output_device,'(a65,/)') & 
         ' Periodic boundary conditions: the minimal image technics is used'

  end subroutine write_cell_vect_to_output
  !******************************************************************

  !******************************************************************
  subroutine cell2vect(ac,bc,cc,alphac,betac,gammac,vc1,vc2,vc3)

    real(kind=r8_kind), intent(in) :: ac,bc,cc,alphac,betac,gammac
    real(kind=r8_kind), intent(out) :: vc1(3),vc2(3),vc3(3)

    real(kind=r8_kind) :: alpha_r, beta_r, gamma_r

    alpha_r=pi*alphac/pi_degree
    beta_r=pi*betac/pi_degree
    gamma_r=pi*gammac/pi_degree

    vc1(1)=ac; vc1(2)=zero; vc1(3)=zero
    vc2(1)=bc*cos(gamma_r); vc2(2)=bc*sin(gamma_r); vc2(3)=zero
    vc3(1)=cc*cos(beta_r)
    vc3(2)=(bc*cc*cos(alpha_r)-vc2(1)*vc3(1))/vc2(2)
    vc3(3)=sqrt(cc*cc-vc3(1)*vc3(1)-vc3(2)*vc3(2))

  end subroutine cell2vect
  !******************************************************************

  !******************************************************************
  subroutine vect2cell(vc1,vc2,vc3,ac,bc,cc,alphac,betac,gammac)

    real(kind=r8_kind), intent(in) :: vc1(3),vc2(3),vc3(3)
    real(kind=r8_kind), intent(out) :: ac,bc,cc,alphac,betac,gammac

    real(kind=r8_kind) :: alpha_r, beta_r, gamma_r
    integer(kind=i4_kind) :: i

    ac=sqrt(dot_product(vc1,vc1))
    bc=sqrt(dot_product(vc2,vc2))
    cc=sqrt(dot_product(vc3,vc3))

    alpha_r=zero; beta_r=zero; gamma_r=zero
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

  !******************************************************************
  function cell_volume()

    real(r8_kind) :: cell_volume
    real(r8_kind) :: v(3,3)
 
      v(1,:)=vect%v1
      v(2,:)=vect%v2
      v(3,:)=vect%v3

      cell_volume=abs(v(1,1)*(v(2,2)*v(3,3)-v(2,3)*v(3,2))+ &
           v(1,2)*(v(3,1)*v(2,3)-v(3,3)*v(2,1))+ &
           v(1,3)*(v(2,1)*v(3,2)-v(2,2)*v(3,1)))

  end function cell_volume
  !******************************************************************

  !******************************************************************
  subroutine calc_k_vectors(vc1,vc2,vc3,v,kv1,kv2,kv3)

    real(kind=r8_kind), intent(in) :: vc1(3),vc2(3),vc3(3),v
    real(kind=r8_kind), intent(out) :: kv1(3),kv2(3),kv3(3)
    
    kv1=two*pi*vector_product(vc2,vc3)/v
    kv2=two*pi*vector_product(vc3,vc1)/v
    kv3=two*pi*vector_product(vc1,vc2)/v

  end subroutine calc_k_vectors
  !******************************************************************

  !******************************************************************
  function valid_cutoff(rcut)

    real(r8_kind) :: valid_cutoff
    real(r8_kind),intent(in) :: rcut
    real(r8_kind) :: rmax
    real(r8_kind) :: amax,bmax,cmax,a_sin,b_sin,g_sin

    a_sin=sin(pi*cel%alpha/pi_degree)
    b_sin=sin(pi*cel%beta/pi_degree)
    g_sin=sin(pi*cel%gamma/pi_degree)
    
    amax=cel%a*b_sin*g_sin
    bmax=cel%b*g_sin*a_sin
    cmax=cel%c*a_sin*b_sin

    rmax=half*min(amax,bmax,cmax)

    if(rcut <= rmax) then
       valid_cutoff=zero
    else
       valid_cutoff=rmax
    end if

  end function valid_cutoff
  !******************************************************************

  !******************************************************************
  subroutine image(r)
    
    real(r8_kind), intent(inout) :: r(3)

    real(r8_kind) :: fr(3),bmat(3,3),det,rr
    integer(i4_kind) :: i

    bmat(1,1)=vect%v2(2)*vect%v3(3)-vect%v2(3)*vect%v3(2)
    bmat(2,1)=vect%v1(3)*vect%v3(2)-vect%v1(2)*vect%v3(3)
    bmat(3,1)=vect%v1(2)*vect%v2(3)-vect%v1(3)*vect%v2(2)
    bmat(1,2)=vect%v2(3)*vect%v3(1)-vect%v2(1)*vect%v3(3)
    bmat(2,2)=vect%v1(1)*vect%v3(3)-vect%v1(3)*vect%v3(1)
    bmat(3,2)=vect%v1(3)*vect%v2(1)-vect%v1(1)*vect%v2(3)
    bmat(1,3)=vect%v2(1)*vect%v3(2)-vect%v2(2)*vect%v3(1)
    bmat(2,3)=vect%v1(2)*vect%v3(1)-vect%v1(1)*vect%v3(2)
    bmat(3,3)=vect%v1(1)*vect%v2(2)-vect%v1(2)*vect%v2(1)

    det=vect%v1(1)*bmat(1,1)+vect%v2(1)*bmat(2,1)+vect%v3(1)*bmat(3,1)
    rr=zero
    if(abs(det) > zero) rr=one/det

    bmat=bmat*rr

    fr=matmul(bmat,r)
    do i=1,3
       if(abs(fr(i)) <= half-small1 .or. abs(fr(i)) >= half+small1) then 
          fr(i)=fr(i)-anint(fr(i))
       end if
    end do
    
    bmat(:,1)=vect%v1
    bmat(:,2)=vect%v2
    bmat(:,3)=vect%v3

    r=matmul(bmat,fr)
   
  end subroutine image
  !******************************************************************

  !******************************************************************
  subroutine calc_images_species(rcut)
    !calculate coordinates of species images
    !not only into nearest cells
    real(r8_kind),intent(in) :: rcut
    real(r8_kind) :: amax,bmax,cmax,a_sin,b_sin,g_sin
    integer(i4_kind) :: na,nb,nc
    integer(i4_kind) :: status,i,j,k,l,m,n
    integer(i4_kind) :: n_image_coor
    real(r8_kind) :: bmat(3,3)

    a_sin=sin(pi*cel%alpha/pi_degree)
    b_sin=sin(pi*cel%beta/pi_degree)
    g_sin=sin(pi*cel%gamma/pi_degree)
    
    amax=cel%a*b_sin*g_sin*half
    bmax=cel%b*g_sin*a_sin*half
    cmax=cel%c*a_sin*b_sin*half

    na=ceiling(rcut/amax); nb=ceiling(rcut/bmax); nc=ceiling(rcut/cmax)
    if(mod(na,2) == 0) then 
       na=na+1
    else
       na=na+2
    end if
    if(mod(nb,2) == 0) then 
       nb=nb+1
    else
       nb=nb+2
    end if
    if(mod(nc,2) == 0) then 
       nc=nc+1
    else
       nc=nc+2
    end if

    n_images=na*nb*nc-1
    n_image_coor=n_images*n_species
    if(.not.allocated(im_coor)) then
!!$       allocate(im_coor(n_image_coor),stat=status)
       allocate(im_coor(n_species),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR allocation")
       do i=1,n_species
          nullify(im_coor(i)%r)
       end do
    end if

    na=(na-1)/2; nb=(nb-1)/2; nc=(nc-1)/2

    m=0
    l_ns: do l=1,n_species
       if(associated(im_coor(l)%r)) then
          deallocate(im_coor(l)%r,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed IM_COOR%R deallocation")
          nullify(im_coor(l)%r)
       end if
       allocate(im_coor(l)%r(3,n_images),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR%R allocation")
       n=0
       l_na: do i=-na,na
          l_nb: do j=-nb,nb
             l_nc: do k=-nc,nc
                if(i == 0 .and. j == 0 .and. k == 0) cycle l_nc
                m=m+1
                n=n+1
                im_coor(l)%r(1,n)=atoms_frac(l)%r(1)+i
                im_coor(l)%r(2,n)=atoms_frac(l)%r(2)+j
                im_coor(l)%r(3,n)=atoms_frac(l)%r(3)+k
                im_coor(l)%type=atoms_cart(l)%type
             end do l_nc
          end do l_nb
       end do l_na
    end do l_ns
    if(m /= n_image_coor) call error_handler( &
         "MolMech: something wrong with image atom")

    !converting fractional coordinates to cartesian
    bmat(:,1)=vect%v1
    bmat(:,2)=vect%v2
    bmat(:,3)=vect%v3

    do j=1,n_species
       do i=1,n_images
          im_coor(j)%r(:,i)=matmul(bmat,im_coor(j)%r(:,i))
       end do
    end do

  end subroutine calc_images_species
  !******************************************************************

  !******************************************************************
  subroutine calc_images_species_solvation(rcut)
    !calculate coordinates of species images in nearest 27 (3x3x3-1) cells
    real(r8_kind),intent(in) :: rcut
    real(r8_kind) :: amax,bmax,cmax,a_sin,b_sin,g_sin
    integer(i4_kind) :: na,nb,nc
    integer(i4_kind) :: status,i,j,k,l,m,n
    integer(i4_kind) :: n_image_coor
    real(r8_kind) :: bmat(3,3)

    a_sin=sin(pi*cel%alpha/pi_degree)
    b_sin=sin(pi*cel%beta/pi_degree)
    g_sin=sin(pi*cel%gamma/pi_degree)
    
    amax=cel%a*b_sin*g_sin*half
    bmax=cel%b*g_sin*a_sin*half
    cmax=cel%c*a_sin*b_sin*half

    na=ceiling(rcut/amax); nb=ceiling(rcut/bmax); nc=ceiling(rcut/cmax)
    if(mod(na,2) == 0) then 
       na=na+1
    else
       na=na+2
    end if
    if(mod(nb,2) == 0) then 
       nb=nb+1
    else
       nb=nb+2
    end if
    if(mod(nc,2) == 0) then 
       nc=nc+1
    else
       nc=nc+2
    end if

    n_images_solv=na*nb*nc-1
    n_image_coor=n_images_solv*n_species
    if(.not.allocated(im_coor_solv)) then
       allocate(im_coor_solv(n_species),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR_SOLV 3 allocation")
       do i=1,n_species
          nullify(im_coor_solv(i)%r)
       end do
    end if

    na=(na-1)/2; nb=(nb-1)/2; nc=(nc-1)/2

    m=0
    l_ns: do l=1,n_species
       if(associated(im_coor_solv(l)%r)) then
          deallocate(im_coor_solv(l)%r,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed IM_COOR_SOLV%R 3 deallocation")
          nullify(im_coor_solv(l)%r)
       end if
       allocate(im_coor_solv(l)%r(3,n_images_solv),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR_SOLV%R 3 allocation")
       n=0
       l_na: do i=-na,na
          l_nb: do j=-nb,nb
             l_nc: do k=-nc,nc
                if(i == 0 .and. j == 0 .and. k == 0) cycle l_nc
                m=m+1
                n=n+1
                im_coor_solv(l)%r(1,n)=atoms_frac(l)%r(1)+i
                im_coor_solv(l)%r(2,n)=atoms_frac(l)%r(2)+j
                im_coor_solv(l)%r(3,n)=atoms_frac(l)%r(3)+k
                im_coor_solv(l)%type=atoms_cart(l)%type
             end do l_nc
          end do l_nb
       end do l_na
    end do l_ns
    ASSERT(m == n_image_coor)

    !converting fractional coordinates to cartesian
    bmat(:,1)=vect%v1
    bmat(:,2)=vect%v2
    bmat(:,3)=vect%v3

    do j=1,n_species
       do i=1,n_images_solv
          im_coor_solv(j)%r(:,i)=matmul(bmat,im_coor_solv(j)%r(:,i))
       end do
    end do

  end subroutine calc_images_species_solvation
  !******************************************************************

  !******************************************************************
  subroutine renew_lattice(new_strain)

    real(r8_kind) :: new_strain(6)
    real(r8_kind) :: strain_mat(3,3)
    real(r8_kind) :: v1(3),v2(3),v3(3),v1p(3),v2p(3),v3p(3)
    integer(i4_kind) :: i

    strain_mat(1,1)=one+new_strain(1)
    strain_mat(2,2)=one+new_strain(2)
    strain_mat(3,3)=one+new_strain(3)
    strain_mat(1,2)=half*new_strain(6)
    strain_mat(2,1)=strain_mat(1,2)
    strain_mat(1,3)=half*new_strain(5)
    strain_mat(3,1)=strain_mat(1,3)
    strain_mat(2,3)=half*new_strain(4)
    strain_mat(3,2)=strain_mat(2,3)

    v1p=vect%v1; v2p=vect%v2; v3p=vect%v3

    v1=zero; v2=zero; v3=zero
    do i=1,3
       v1(i)=v1(i)+strain_mat(1,i)*v1p(1)+strain_mat(2,i)*v1p(2)+strain_mat(3,i)*v1p(3)
       v2(i)=v2(i)+strain_mat(1,i)*v2p(1)+strain_mat(2,i)*v2p(2)+strain_mat(3,i)*v2p(3)
       v3(i)=v3(i)+strain_mat(1,i)*v3p(1)+strain_mat(2,i)*v3p(2)+strain_mat(3,i)*v3p(3)
    end do

    vect%v1=v1; vect%v2=v2; vect%v3=v3

    call vect2cell(vect%v1,vect%v2,vect%v3, &
         cel%a,cel%b,cel%c,cel%alpha,cel%beta,cel%gamma)
    volume=cell_volume()
    call calc_k_vectors(vect%v1,vect%v2,vect%v3,volume,kvect%v1,kvect%v2,kvect%v3)

  end subroutine renew_lattice
  !******************************************************************

  !******************************************************************
  subroutine store_vect()

    vect_store=vect

  end subroutine store_vect
  !******************************************************************

  !******************************************************************
  subroutine use_vect()

    vect=vect_store

  end subroutine use_vect
  !******************************************************************

  !******************************************************************
  function vector_product(v1,v2)

    real(kind=r8_kind) :: vector_product(3)
    real(kind=r8_kind) :: v1(3),v2(3)

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product
  !******************************************************************

  !******************************************************************
  subroutine calc_images_slab_species(rcut)
    !calculate coordinates of species images
    !not only into nearest slab cells
    real(r8_kind),intent(in) :: rcut
    real(r8_kind) :: amax,bmax,a_sin
    integer(i4_kind) :: na,nb
    integer(i4_kind) :: status,i,j,l,m,n
    integer(i4_kind) :: n_image_coor
    real(r8_kind) :: bmat(2,2)

    a_sin=sin(pi*cel_s%alpha/pi_degree)
    
    amax=cel_s%a*a_sin*half
    bmax=cel_s%b*a_sin*half

    na=ceiling(rcut/amax); nb=ceiling(rcut/bmax)
    if(mod(na,2) == 0) then 
       na=na+1
    else
       na=na+2
    end if
    if(mod(nb,2) == 0) then 
       nb=nb+1
    else
       nb=nb+2
    end if

    n_images=na*nb-1
    n_image_coor=n_images*n_species
    if(.not.allocated(im_coor)) then
!!$       allocate(im_coor(n_image_coor),stat=status)
       allocate(im_coor(n_species),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR allocation")
       do i=1,n_species
          nullify(im_coor(i)%r)
       end do
    end if

    na=(na-1)/2; nb=(nb-1)/2

    m=0
    l_ns: do l=1,n_species
       if(associated(im_coor(l)%r)) then
          deallocate(im_coor(l)%r,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed IM_COOR%R deallocation")
          nullify(im_coor(l)%r)
       end if
       allocate(im_coor(l)%r(3,n_images),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR%R allocation")
       n=0
       l_na: do i=-na,na
          l_nb: do j=-nb,nb
             if(i == 0 .and. j == 0) cycle l_nb
             m=m+1
             n=n+1
             im_coor(l)%r(1,n)=atoms_frac(l)%r(1)+i
             im_coor(l)%r(2,n)=atoms_frac(l)%r(2)+j
             im_coor(l)%r(3,n)=atoms_frac(l)%r(3)
             im_coor(l)%type=atoms_cart(l)%type
          end do l_nb
       end do l_na
    end do l_ns
    if(m /= n_image_coor) call error_handler( &
         "MolMech: something wrong with image atom")

    !converting fractional coordinates to cartesian
    bmat(:,1)=vect_s%v1
    bmat(:,2)=vect_s%v2

    do j=1,n_species
       do i=1,n_images
          im_coor(j)%r(1:2,i)=matmul(bmat,im_coor(j)%r(1:2,i))
       end do
    end do

  end subroutine calc_images_slab_species
  !******************************************************************

  !******************************************************************
  subroutine calc_images_slab_species_solvation(rcut)
    !calculate coordinates of species images in nearest 8 (3x3-1) cells
    real(r8_kind),intent(in) :: rcut

    real(r8_kind) :: amax,bmax,a_sin
    integer(i4_kind) :: na,nb
    integer(i4_kind) :: status,i,j,l,m,n
    integer(i4_kind) :: n_image_coor
    real(r8_kind) :: bmat(2,2)

    a_sin=sin(pi*cel_s%alpha/pi_degree)
    
    amax=cel_s%a*a_sin*half
    bmax=cel_s%b*a_sin*half

    na=ceiling(rcut/amax); nb=ceiling(rcut/bmax)
    if(mod(na,2) == 0) then 
       na=na+1
    else
       na=na+2
    end if
    if(mod(nb,2) == 0) then 
       nb=nb+1
    else
       nb=nb+2
    end if

!    na=3; nb=3

    n_images_solv=na*nb-1
    n_image_coor=n_images_solv*n_species
    if(.not.allocated(im_coor_solv)) then
       allocate(im_coor_solv(n_species),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR_SOLV allocation")
       do i=1,n_species
          nullify(im_coor_solv(i)%r)
       end do
    end if

    na=(na-1)/2; nb=(nb-1)/2

    m=0
    l_ns: do l=1,n_species
       if(associated(im_coor_solv(l)%r)) then
          deallocate(im_coor_solv(l)%r,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed IM_COOR%R deallocation")
          nullify(im_coor_solv(l)%r)
       end if
       allocate(im_coor_solv(l)%r(3,n_images_solv),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR_SOLV%R allocation")
       n=0
       l_na: do i=-na,na
          l_nb: do j=-nb,nb
             if(i == 0 .and. j == 0) cycle l_nb
             m=m+1
             n=n+1
             im_coor_solv(l)%r(1,n)=atoms_frac(l)%r(1)+i
             im_coor_solv(l)%r(2,n)=atoms_frac(l)%r(2)+j
             im_coor_solv(l)%r(3,n)=atoms_frac(l)%r(3)
             im_coor_solv(l)%type=atoms_cart(l)%type
          end do l_nb
       end do l_na
    end do l_ns
    ASSERT(m == n_image_coor)

    !converting fractional coordinates to cartesian
    bmat(:,1)=vect_s%v1
    bmat(:,2)=vect_s%v2

    do j=1,n_species
       do i=1,n_images_solv
          im_coor_solv(j)%r(1:2,i)=matmul(bmat,im_coor_solv(j)%r(1:2,i))
       end do
    end do

  end subroutine calc_images_slab_species_solvation
  !******************************************************************

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !******************************************************************
  function read_coord()

    logical ::read_coord
    integer(kind=i4_kind) :: status,i,j,i0,i1,j0,j1,ls
    logical :: cart_s,cart_f,coor_cart, frac_s,frac_f,coor_frac
    character(len=16) :: coor_type
    character(len=80) :: string
    character(len=6) :: number
    character(len=len_name) :: a_nm
    character(len=1) :: c_s
    real(kind=r8_kind) :: coor(3)
    integer(kind=i4_kind), allocatable :: at_type(:), at_cs(:)
    character(len=len_name), allocatable :: at_nm(:)

    n_fixed=0

    if(calc_strain) then
       if(lattice_calc) strain=one
       if(slab_calc) strain_s=one
    end if

    if(lattice_calc) then
       ls=3
    elseif(slab_calc) then
       ls=2
    end if

    !Are coordinates cartesian?
    cart_s=find_word("&CARTESIAN_COOR",i0)
    cart_f=find_word("/CARTESIAN_COOR",i1)
    coor_cart=(cart_s .and. cart_f)
    if(.not. coor_cart) then
       if(.not.cart_s .and. cart_f) call error_handler &
            ("MolMech: You specified Cartesian coordinates block without &CARTESIAN_COOR")
       if(cart_s .and. .not.cart_f) call error_handler &
            ("MolMech: You specified Cartesian coordinates block without /CARTESIAN_COOR")
    end if

    !Are coordinates fractional?
    frac_s=find_word("&FRACTIONAL_COOR",j0)
    frac_f=find_word("/FRACTIONAL_COOR",j1)
    coor_frac=(frac_s .and. frac_f)
    if(.not. coor_frac) then
       if(.not.frac_s .and. frac_f) call error_handler &
            ("MolMech: You specified Fractional coordinates block without &FRACTIONAL_COOR")
       if(frac_s .and. .not.frac_f) call error_handler &
            ("MolMech: You specified Fractional coordinates block without /FRACTIONAL_COOR")
    end if

    if(.not. coor_cart .and. .not. coor_frac) then
       read_coord=.false.
       return
    end if

    if(coor_cart .and. coor_frac) call error_handler &
         ("MolMech: Pleace specify either fractional or cartesian coordinates ")
    
    if(coor_cart) then
       n_species=i1-i0-1
       coor_type="&CARTESIAN_COOR"
    end if
    if(coor_frac) then
       n_species=j1-j0-1
       coor_type="&FRACTIONAL_COOR"
    end if

    if(n_species<=0)  then
       read_coord=.false.
       return
    end if
    
    allocate(atoms_cart(n_species),stat=status)
    if(status /= 0) call error_handler( &
         "MolMech: failed ATOMS_CART allocation")

    if(lattice_calc .or. slab_calc) then
       allocate(atoms_frac(n_species),stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed ATOMS_FRAC allocation")
    end if

    if(coor_frac .and. .not.(lattice_calc .or. slab_calc)) then
       call error_handler( &
            "MolMech: You are going to perform non periodic calculation. &
           & Use cartesian coordinates!")
    end if

    allocate(at_type(n_species),at_cs(n_species),at_nm(n_species), stat=status)
    if(status /= 0) call error_handler( &
            "MolMech: failed AT_TYPE allocation")

    n_species_types=0_i4_kind

    read_coord=find_word(trim(coor_type),i)
    do i=1,n_species
       read(input_device,'(a80)', err=300) string
       if(check_string(string," c ") .or. &
          check_string(string," C ") .or. &
          check_string(string," s ") .or. &
          check_string(string," S ")) then
          read(string, *, err=300) a_nm,c_s,coor(1:3)
          call upcase(c_s)
       else
          read(string, *,err=300) a_nm,coor(1:3)
          c_s="C"
       end if
       if(coor_cart) then
          atoms_cart(i)%type=get_atom_type()
          atoms_cart(i)%initial_number=i
          atoms_cart(i)%neigh_sp=i
          atoms_cart(i)%r=coor
       else 
          atoms_frac(i)%type=get_atom_type()
          atoms_frac(i)%initial_number=i
          atoms_frac(i)%r=coor
          do j=1,ls
             if(atoms_frac(i)%r(j) < zero) then 
                atoms_frac(i)%r(j)=one+atoms_frac(i)%r(j)
             elseif(atoms_frac(i)%r(j) >= one) then
                atoms_frac(i)%r(j)=atoms_frac(i)%r(j)-one
             end if
          end do
          atoms_frac(i)%neigh_sp=i
       end if
    end do

    if(lattice_calc) then
       select case (coor_cart)
          case (.true.)
             call cart2frac(n_species,vect,atoms_cart,atoms_frac)
             call frac2cart(n_species,vect,atoms_frac,atoms_cart)
          case (.false.)
             call frac2cart(n_species,vect,atoms_frac,atoms_cart)
             do i=1,n_species
                atoms_cart(i)%type=atoms_frac(i)%type
                atoms_cart(i)%initial_number=atoms_frac(i)%initial_number
                atoms_cart(i)%neigh_sp=atoms_frac(i)%neigh_sp
             end do
       end select
    else if(slab_calc) then
       select case (coor_cart)
          case (.true.)
             do i=1,n_species
                call cart2frac_slab(atoms_cart(i)%r,atoms_frac(i)%r)
                call frac2cart_slab(atoms_frac(i)%r,atoms_cart(i)%r)
             end do
          case (.false.)
             do i=1,n_species
                call frac2cart_slab(atoms_frac(i)%r,atoms_cart(i)%r)
                atoms_cart(i)%type=atoms_frac(i)%type
                atoms_cart(i)%initial_number=atoms_frac(i)%initial_number
                atoms_cart(i)%neigh_sp=atoms_frac(i)%neigh_sp
             end do
       end select
    end if

    allocate(atoms(n_species_types), stat=status)
    if(status /= 0) call error_handler( &
            "MolMech: failed ATOMS allocation")

    do i=1,n_species_types
       atoms(i)%name=at_nm(i)
       atoms(i)%c_s=at_cs(i)
    end do

    deallocate(at_type,at_cs,at_nm, stat=status)
    if(status /= 0) call error_handler( &
            "MolMech: failed AT_TYPE deallocation")
    read_coord=.true.
    return

300 write(number,'(i4)') i
    call error_handler("MolMech: Coordinates of "//trim(number)// &
         " atom has been read in with error")

  contains
    
    function get_atom_type()

      integer(kind=i4_kind) :: get_atom_type

      integer(kind=i4_kind) :: ig,csg

      if(n_species_types==0) then
         n_species_types=1
         at_nm(1)=a_nm
         at_type(1)=n_species_types
         if(c_s == "C") at_cs(1)=0
         if(c_s == "S") at_cs(1)=1
         get_atom_type=at_type(1)
      else
         do ig=1,n_species_types
            if(c_s == "C") csg=0
            if(c_s == "S") csg=1
            if(a_nm == at_nm(ig) .and. csg == at_cs(ig)) then
               get_atom_type=at_type(ig)
               return
            end if
         end do
         n_species_types=n_species_types+1
         at_nm(n_species_types)=a_nm
         at_type(n_species_types)=n_species_types
         if(c_s == "C") at_cs(n_species_types)=0
         if(c_s == "S") at_cs(n_species_types)=1
         get_atom_type=at_type(n_species_types)
      end if

    end function get_atom_type

  end function   read_coord
  !******************************************************************

  !******************************************************************
  function read_coord_gx()

    logical ::read_coord_gx

    integer :: gx_file
    real(r8_kind) :: atnum, coor(3)
    integer :: unum, num, znum(3), cnum(3)
    integer :: i
    character(len=6) :: number

    if(with_optimizer) then
       call get_file_device(gx_file,'gxfile','inp')
    else
       call get_file_device(gx_file,'gx.mm','inp')
    end if

    i=0
    loop_i: do 
       if(i == n_species) exit loop_i
       read(gx_file,*,end=100,err=110) atnum,coor,unum,num,znum,cnum
       if(atnum == gx_dummy_atom) cycle loop_i
       if(atnum < zero .and. i < n_species) goto 100
       i=i+1
       atoms_cart(i)%r = coor*b2a
    end do loop_i
    read(gx_file,*,err=110) gx_work

    if(lattice_calc) then
       call cart2frac(n_species,vect,atoms_cart,atoms_frac)
       call frac2cart(n_species,vect,atoms_frac,atoms_cart)
    else if(slab_calc) then 
       do i=1,n_species
          call cart2frac_slab(atoms_cart(i)%r,atoms_frac(i)%r)
          call frac2cart_slab(atoms_frac(i)%r,atoms_cart(i)%r)
       end do
    end if

    read_coord_gx = .true.
    call close_file_device(gx_file)
    return

100 read_coord_gx = .false.
    call close_file_device(gx_file)
    return

110 write(number,'(i4)') i
    call error_handler("MolMech: Coordinates of "//trim(number)// &
         "atom has been read in from GX.MM file with error")

  end function read_coord_gx
  !******************************************************************

  !******************************************************************
  function read_qmmm()
    !for IMOMM

    logical :: read_qmmm
    integer(i4_kind) :: i,j

    if(qmmm == 2) then
       j=0
       do i=1,n_qm_at
          if(gx_qm(i)%charge == 99.0_r8_kind) cycle
          j=j+1
          atoms_cart(j)%r(1)=gx_qm(i)%x*b2a
          atoms_cart(j)%r(2)=gx_qm(i)%y*b2a
          atoms_cart(j)%r(3)=gx_qm(i)%z*b2a
       end do
       if(n_species /= j) call error_handler( &
            "MolMech: Number of atoms in molmech.inp.2 does not coincide with "// &
            "calculated by qmmm interface")
    elseif(qmmm == 3) then
       j=0
       do i=1,n_total_at
          if(gx(i)%charge == 99.0_r8_kind) cycle
          j=j+1
          atoms_cart(j)%r(1)=gx(i)%x*b2a
          atoms_cart(j)%r(2)=gx(i)%y*b2a
          atoms_cart(j)%r(3)=gx(i)%z*b2a
       end do
       if(n_species /= j) call error_handler( &
            "MolMech: Number of atoms in molmech.inp.3 does not coincide with "// &
            "calculated by qmmm interface")
    end if

    if(lattice_calc) then
       call cart2frac(n_species,vect,atoms_cart,atoms_frac)
       call frac2cart(n_species,vect,atoms_frac,atoms_cart)
    else if(slab_calc) then 
       do i=1,n_species
          call cart2frac_slab(atoms_cart(i)%r,atoms_frac(i)%r)
          call frac2cart_slab(atoms_frac(i)%r,atoms_cart(i)%r)
       end do
    end if

    read_qmmm = .true.

  end function read_qmmm
  !******************************************************************

  !******************************************************************
  subroutine cart2frac(natms,amat,cart,frac)

    integer(kind=i4_kind), intent(in) :: natms
    type(lat_vectors) ,intent(in) :: amat
    type(atom_coord), intent(in) :: cart(natms)
    type(atom_coord), intent(out) :: frac(natms)

    real(kind=r8_kind) :: bmat(3,3),det,rr
    integer(kind=i4_kind) :: i,j

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
    rr=zero
    if(abs(det) > zero) rr=one/det

    bmat=bmat*rr

    do i=1,natms
       frac(i)%r=matmul(bmat,cart(i)%r)
       do j=1,3
          if(abs(frac(i)%r(j)) <= half-small1 .or. abs(frac(i)%r(j)) >= half+small1) then
             frac(i)%r(j)=frac(i)%r(j)-anint(frac(i)%r(j))
          end if
          if(abs(frac(i)%r(j)) <= small1) frac(i)%r(j)=zero
          if(frac(i)%r(j) < zero) then 
             frac(i)%r(j)=one+frac(i)%r(j)
          elseif(frac(i)%r(j) >= one) then
             frac(i)%r(j)=frac(i)%r(j)-one
          end if
       end do
    end do

  end subroutine cart2frac
  !******************************************************************

  !******************************************************************
  subroutine frac2cart(natms,amat,frac,cart)

    integer(kind=i4_kind), intent(in) :: natms
    type(lat_vectors) ,intent(in) :: amat
    type(atom_coord), intent(in) :: frac(natms)
    type(atom_coord)  :: cart(natms)

    real(kind=r8_kind) :: bmat(3,3)
    integer(kind=i4_kind) :: i

    bmat(:,1)=amat%v1
    bmat(:,2)=amat%v2
    bmat(:,3)=amat%v3

    do i=1,natms
       cart(i)%r=matmul(bmat,frac(i)%r)
    end do

  end subroutine frac2cart
  !******************************************************************

  !******************************************************************
  subroutine read_species()

    integer(i4_kind) :: i,j,k,status,ii
    character(len=6) :: nm
    character(len=4) :: number
    character(len=1) :: c_s
    character(len=2) :: buf
    logical, allocatable :: exist(:)
    logical :: read_spec

    atoms(1:n_species_types)%charge=df_charge
    atoms(1:n_species_types)%redfac=df_vdw_red_fac
    atoms(1:n_species_types)%r_coval=df_r_coval
    
    allocate(exist(n_species_types), stat=status)
    if(status /= 0) call error_handler( &
            "MolMech: failed EXIST allocation")
    exist=.false.

    call go_to_first_input_line
!!$    lab1: do i=1,n_species_types
    lab1: do 
       charge=df_charge
       r_coval=df_r_coval
       r_vdw=df_r_vdw
       r0_dr=df_r0_dr
       k_dr=df_k_dr
       name=df_name
       main_name=df_main_name
       vdw_red_fac=df_vdw_red_fac
       solv_scale_fac=df_solv_scale_fac
       solv_scale_fac_smooth=df_solv_scale_fac_smooth
       no_sphere=df_no_sphere

       read_spec=find_namelist("&SPECIES",ii)
       if(.not.read_spec) goto 100
       read(input_device,nml=species, end=100, err=200)
       if(trim(name)==" ") goto 300
       if(check_string(name," s").or.check_string(name," S")) then
          c_s="S"
          call name_without_cs(name,nm)
       elseif(check_string(name," c").or.check_string(name," C")) then
          c_s="C"
          call name_without_cs(name,nm)
       else
          c_s="C"
          nm=name
       end if
       lab2: do j=1,n_species_types
          if(j==name2type(nm,c_s)) then
             if(exist(j)) goto 400
             exist(j)=.true.
             atoms(j)%charge=charge
             atoms(j)%redfac=vdw_red_fac
             if(r_coval==df_r_coval) then
                buf=nm(1:2)
                call upcase(buf)
                call upcase(main_name)
                lab3: do k=1,n_elements
                   if(atom_data(k)%name == buf .or. &
                        atom_data(k)%name == main_name) then
                      r_coval=atom_data(k)%r_coval
                      exit
                   end if
                end do lab3
             end if
             atoms(j)%r_coval=r_coval
             if(r_vdw==df_r_vdw) then
                buf=nm(1:2)
                call upcase(buf)
                call upcase(main_name)
                lab4: do k=1,n_elements
                   if(atom_data(k)%name == buf .or. &
                        atom_data(k)%name == main_name) then
                      r_vdw=atom_data(k)%r_vdw_bondi
                      exit
                   end if
                end do lab4
             end if
             atoms(j)%r_vdw=r_vdw
             if(r0_dr==df_r0_dr) then
                buf=nm(1:2)
                call upcase(buf)
                call upcase(main_name)
                lab5: do k=1,n_elements
                   if(atom_data(k)%name == buf .or. &
                        atom_data(k)%name == main_name) then
                      r0_dr=atom_data(k)%r0_dr
                      exit
                   end if
                end do lab5
             end if
             atoms(j)%r0_dr=r0_dr
             if(k_dr==df_k_dr) then
                buf=nm(1:2)
                call upcase(buf)
                call upcase(main_name)
                lab6: do k=1,n_elements
                   if(atom_data(k)%name == buf .or. &
                        atom_data(k)%name == main_name) then
                      k_dr=atom_data(k)%k_dr
                      exit
                   end if
                end do lab6
             end if
             atoms(j)%k_dr=k_dr
             if(solv_scale_fac==df_solv_scale_fac) then
                buf=nm(1:2)
                call upcase(buf)
                call upcase(main_name)
                lab7: do k=1,n_elements
                   if(atom_data(k)%name == buf .or. &
                        atom_data(k)%name == main_name) then
                      solv_scale_fac=atom_data(k)%solv_scale_fac
                      exit
                   end if
                end do lab7
             end if
             atoms(j)%scalefac=solv_scale_fac
             if(solv_scale_fac_smooth==df_solv_scale_fac_smooth) then
                buf=nm(1:2)
                call upcase(buf)
                call upcase(main_name)
                lab8: do k=1,n_elements
                   if(atom_data(k)%name == buf .or. &
                        atom_data(k)%name == main_name) then
                      solv_scale_fac_smooth=atom_data(k)%solv_scale_fac_smooth
                      exit
                   end if
                end do lab8
             end if
             atoms(j)%scalefac_smooth=solv_scale_fac_smooth
             atoms(j)%no_sphere=no_sphere
             if(atoms(j)%c_s == 1) atoms(j)%no_sphere=.true. !no spheres for shells
             cycle lab1
          end if
       end do lab2
!!$       goto 500
    end do lab1

100 deallocate(exist, stat=status)
    if(status /= 0) call error_handler( &
            "MolMech: failed EXIST deallocation")

    return

!!$100 write(number,'(i4)') n_species_types
!!$    call error_handler("MolMech: Your input has to have "// &
!!$         trim(number)//" namelists SPECIES")

200 call input_nm_error(i,"SPECIES")

300 write(number,'(i4)') i
    call error_handler( &
         "MolMech: Species name has to be specified explisitly in "// &
         trim(number)//" namelist SPECIES")

400 write(number,'(i4)') i
    call error_handler( &
         "MolMech: Repeated definition of species name in "// &
         trim(number)//" namelist SPECIES")

500 write(number,'(i4)') i
    call error_handler( &
         "MolMech: Species name in "//trim(number)// &
         " namelist SPECIES no corresponds to any atom")

  end subroutine read_species
  !******************************************************************

  !******************************************************************
  subroutine read_pdc()

    integer :: pdc_file
    integer :: i,j,typ
    real(r8_kind) :: pdc

    call get_file_device(pdc_file,'pdc.save','inp')

    i=0
    do 
       n_fixed=i
       read(pdc_file,*,end=100) j,pdc
       i=i+1
       if(i > n_species) exit
       typ=atoms_cart(i)%type
       atoms(typ)%charge=pdc
    end do

100 call close_file_device(pdc_file)

  end subroutine read_pdc
  !******************************************************************

  !******************************************************************
  function check_total_charge_and_dipmom()

    real(r8_kind) :: check_total_charge_and_dipmom(4)
    real(r8_kind) :: q,dm(3)
    integer(i4_kind) :: i,a_type

    q=zero; dm=zero
    do i=1,n_species
       a_type=atoms_cart(i)%type
       q=q+atoms(a_type)%charge
       dm=dm+atoms(a_type)%charge*atoms_cart(i)%r
    end do

    check_total_charge_and_dipmom(1)=q
    check_total_charge_and_dipmom(2:4)=dm

    write(output_device,'(/)')
    write(output_device,'(1x,a21,f15.8)') "Total charge:        ", q
    write(output_device,'(1x,a21,3f15.8)') "Total dipole moment: ", dm
    write(output_device,'(/)')

  end function check_total_charge_and_dipmom
  !******************************************************************

  !******************************************************************
  subroutine write_species_to_output()

    integer(kind=i4_kind) :: i,j
    character(len=5) :: c_s

    write(output_device,'(/)')
    write(output_device,'(1x,a40,i6,/)') "Total number of atoms (core and shell): ", &
         n_species

    if(lattice_calc .or. slab_calc) then
       write(output_device,'(1x,a23,/)') "Fractional coordinates:"
       write(output_device,'(80("-"))')
       write(output_device,'(a80)') &
            " Number    Type  Name  Core/Shell         X              Y              Z       "
       write(output_device,'(80("-"))')
       do i=1,n_species
          j=atoms_cart(i)%type
          c_s="core "
          if(atoms(j)%c_s == 1) c_s="shell"
          write(output_device,'(i8,3x,i4,2x,a4,4x,a5,5x,3f15.8)') &
               i,j,atoms(j)%name,c_s,atoms_frac(i)%r(1:3)
       end do
       write(output_device,'(80("-"))')
       write(output_device,'(/)')
    end if

    write(output_device,'(1x,a22,/)') "Cartesian coordinates:"
    write(output_device,'(80("-"))')
    write(output_device,'(a80)') &
         " Number    Type  Name  Core/Shell         X              Y              Z       "
    write(output_device,'(80("-"))')
    do i=1,n_species
       j=atoms_cart(i)%type
       c_s="core "
       if(atoms(j)%c_s == 1) c_s="shell"
       write(output_device,'(i8,3x,i4,2x,a4,4x,a5,5x,3f15.8)') &
            i,j,atoms(j)%name,c_s,atoms_cart(i)%r(1:3)
    end do
    write(output_device,'(80("-"))')
    write(output_device,'(/)')

    write(output_device,'(1x,a17,/)') "Type description:"
    write(output_device,'(80("-"))')
    if(.not.solvent) then
       write(output_device,'(a56)') &
            " Type  Name  Core/Shell       Charge      R_coval  R_vdw"
    else
       write(output_device,'(a80)') &
            " Type  Name  Core/Shell       Charge      R_coval  R_vdw  Scale_fac Scale_fac_smooth"
    end if
    write(output_device,'(80("-"))')
    do i=1,n_species_types
       c_s="core "
       if(atoms(i)%c_s == 1_i4_kind) c_s="shell"
    if(.not.solvent) then
       write(output_device,'(1x,i4,2x,a4,4x,a5,5x,f15.8,2x,f4.2,5x,f4.2)') &
            i,atoms(i)%name,c_s,atoms(i)%charge,atoms(i)%r_coval,atoms(i)%r_vdw
    else
       write(output_device,'(1x,i4,2x,a4,4x,a5,5x,f15.8,2x,f4.2,5x,f4.2,5x,f4.2,7x,f6.4)') &
            i,atoms(i)%name,c_s,atoms(i)%charge,atoms(i)%r_coval,atoms(i)%r_vdw,atoms(i)%scalefac, &
            atoms(i)%scalefac_smooth
    end if
    end do
    write(output_device,'(80("-"))')
    write(output_device,'(/)')

  end subroutine write_species_to_output
  !******************************************************************

  !******************************************************************
  subroutine write_final_geometry()

    integer(i4_kind) :: i,j,k
    character(len=1) :: c_s

    write(output_device,'(/)')
    write(output_device,'(80("*"))')
    write(output_device,'(a23)') " Final atom coordinates"
    if(lattice_calc .or. slab_calc) then
       write(output_device,'(a12)') " fractional:"
       write(output_device,'(80("-"))')
       write(output_device,'(a79)') &
            " Number  Init_num  Type Name c/s         X              Y              Z       "
       write(output_device,'(80("-"))')
       do i=1,n_species
          j=atoms_cart(i)%type
          k=atoms_cart(i)%initial_number
          c_s="c"
          if(atoms(j)%c_s == 1) c_s="s"
          write(output_device,'(i8,1x,i8,2x,i4,1x,a6,1x,a1,1x,2x,3f15.8)') &
               i,k,j,atoms(j)%name,c_s,atoms_frac(i)%r(1:3)
       end do
       write(output_device,'(80("-"))')
    end if

    write(output_device,'(a11)') " cartesian:"
    write(output_device,'(80("-"))')
    write(output_device,'(a79)') &
         " Number  Init_num  Type Name c/s         X              Y              Z       "
    write(output_device,'(80("-"))')
    do i=1,n_species
       j=atoms_cart(i)%type
       k=atoms_cart(i)%initial_number
       c_s="c"
       if(atoms(j)%c_s == 1) c_s="s"
       write(output_device,'(i8,1x,i8,2x,i4,1x,a6,1x,a1,1x,2x,3f15.8)') &
            i,k,j,atoms(j)%name,c_s,atoms_cart(i)%r(1:3)
    end do

  end subroutine write_final_geometry
  !******************************************************************

  !******************************************************************
  subroutine save_vm()
    !for Viewmol

    integer(i4_kind) :: file_vm
    integer(i4_kind) :: i,j,k,l,ich
    character(len=len_name) :: a_name
    real(r8_kind) :: max_z,min_z,z_gap

    call get_file_device(file_vm,'molmech.vm','out')

    write(file_vm,'(a7)')  " $title"
    write(file_vm,'(1x,a30)') job_topic
    if(lattice_calc) then 
       write(file_vm,'(a20)')  " $coord fractional 1"
    else
       write(file_vm,'(a10)')  " $coord 1"
    end if
       
    max_z=zero; min_z=zero
    l1: do l=1,n_species
       l2: do i=1,n_species
          if(atoms_cart(i)%initial_number /= l) cycle l2
          j=atoms_cart(i)%type
          if(get_c_s(atoms(j)%c_s) /= "C") cycle l2
          a_name=atoms(j)%name
          do k=1,len_name
             ich=iachar(a_name(k:k))
             if(ich >= 48 .and. ich <= 57) then
                ich=32
                a_name(k:k)=achar(ich)
             end if
          end do
          if(lattice_calc) then 
             write(file_vm,'(3f16.9,1x,a6)') atoms_frac(i)%r(1:3),a_name
          else
             write(file_vm,'(3f16.9,1x,a6)') atoms_cart(i)%r(1:3),a_name
             if(atoms_cart(i)%r(3) < min_z) min_z=atoms_cart(i)%r(3)
             if(atoms_cart(i)%r(3) > max_z) max_z=atoms_cart(i)%r(3)
          end if
          exit l2
       end do l2
    end do l1

    z_gap=(max_z-min_z)*three
 
    if(lattice_calc) write(file_vm,'(a11,6f10.5)') " $unitcell ", &
         cel%a,cel%b,cel%c,cel%alpha,cel%beta,cel%gamma
    if(slab_calc) write(file_vm,'(a11,6f10.5)') " $unitcell ", &
         cel_s%a,cel_s%b,z_gap,90.0_r8_kind,90.0_r8_kind,cel_s%alpha
    write(file_vm,'(a5)')  " $end"

    call close_file_device(file_vm)

  end subroutine save_vm
  !******************************************************************

  !******************************************************************
  subroutine save_xyz()

    integer(i4_kind) :: file_xyz
    integer(i4_kind) :: i,j,k,l,ich,nn
    character(len=len_name) :: a_name
    real(r8_kind) :: v(3)

    call get_file_device(file_xyz,'molmech.xyz','out')

    nn=0
    do i=1,n_species
       j=atoms_cart(i)%type
       if(get_c_s(atoms(j)%c_s) /= "C") cycle
       nn=nn+1
    end do

    write(file_xyz,'(i6)') nn
    write(file_xyz,*) "This file was saved by ParaGauss FF module"

    nn=0
    l1: do l=1,n_species
       l2: do i=1,n_species
          if(atoms_cart(i)%initial_number /= l) cycle l2
          j=atoms_cart(i)%type
          if(get_c_s(atoms(j)%c_s) /= "C") cycle l2
          a_name=atoms(j)%name
          do k=1,len_name
             ich=iachar(a_name(k:k))
             if(ich >= 48 .and. ich <= 57) then
                ich=32
                a_name(k:k)=achar(ich)
             end if
          end do
          nn=nn+1
          if(to_xmkmol_file) then
             if(nn==1) then
                if(lattice_calc .or. slab_calc) then
                   write(file_xyz,'(a6,3f16.9,a29)') &
                        a_name, atoms_cart(i)%r(1:3)," crystal_origin   0.0 0.0 0.0"
                else
                   write(file_xyz,'(a6,3f16.9)') a_name, atoms_cart(i)%r(1:3)
                end if
             else if(nn==2) then
                if(lattice_calc .or. slab_calc) then
                   write(file_xyz,'(a6,3f16.9,a21)') &
                        a_name, atoms_cart(i)%r(1:3)," crystal_images 1 1 1"
                else
                   write(file_xyz,'(a6,3f16.9)') a_name, atoms_cart(i)%r(1:3)
                end if
             else if(nn==3) then
                v=zero
                if(lattice_calc .or. slab_calc) then
                   if(lattice_calc) v=vect%v1
                   if(slab_calc) v(1:2)=vect_s%v1
                   write(file_xyz,'(a6,3f16.9,a18,3f16.9)') &
                        a_name, atoms_cart(i)%r(1:3)," crystal_vector 1 ",v
                else
                   write(file_xyz,'(a6,3f16.9)') a_name, atoms_cart(i)%r(1:3)
                end if
             else if(nn==4) then
                v=zero
                if(lattice_calc .or. slab_calc) then
                   if(lattice_calc) v=vect%v2
                   if(slab_calc) v(1:2)=vect_s%v2
                   write(file_xyz,'(a6,3f16.9,a18,3f16.9)') &
                        a_name, atoms_cart(i)%r(1:3)," crystal_vector 2 ",v
                else
                   write(file_xyz,'(a6,3f16.9)') a_name, atoms_cart(i)%r(1:3)
                end if
             else if(nn==5) then
                v=zero
                if(lattice_calc .or. slab_calc) then
                   if(lattice_calc) v=vect%v3
                   if(slab_calc) v(3)=40.0_r8_kind
                   write(file_xyz,'(a6,3f16.9,a18,3f16.9)') &
                        a_name, atoms_cart(i)%r(1:3)," crystal_vector 3 ",v
                else
                   write(file_xyz,'(a6,3f16.9)') a_name, atoms_cart(i)%r(1:3)
                end if
             else
                write(file_xyz,'(a6,3f16.9)') a_name, atoms_cart(i)%r(1:3)
             end if
          else
             write(file_xyz,'(a6,3f16.9)') a_name, atoms_cart(i)%r(1:3)
          end if
          exit l2
       end do l2
    end do l1

    call close_file_device(file_xyz)

  end subroutine save_xyz
  !******************************************************************

  !******************************************************************
  subroutine epe_interface()

    integer(i4_kind) :: file_gxcv
    integer(i4_kind) :: i,i1,j,k,l,m,ty,lenght,status
    real(r8_kind) :: at,x,y,z,x1,y1,z1
    integer(i4_kind) :: atty(100)
    real(r8_kind), allocatable :: st(:)
    integer (i4_kind), allocatable :: ishell(:)

    call get_file_device(file_gxcv,'cellvec','inp')

    lenght=0
    if(allocated(cs_pair)) then
       lenght=size(cs_pair,2)
       allocate (ishell(lenght), st(lenght), stat=status)
       if(status /= 0) call error_handler("MolMech: failed ISHELL allocation")
    end if

    if(lattice_calc) then
       write(file_gxcv,'(3f15.8)') vect%v1
       write(file_gxcv,'(3f15.8)') vect%v2
       write(file_gxcv,'(3f15.8)') vect%v3
    else if(slab_calc) then
       write(file_gxcv,'(3f15.8)') vect_s%v1,zero
       write(file_gxcv,'(3f15.8)') vect_s%v2,zero
       write(file_gxcv,'(3f15.8)') zero,zero,hundred
    end if

    atty=0
    i1=0
    k=0
    l1: do i=1,N_species
       l2: do j=1,N_species
          if(atoms_cart(j)%initial_number /= i) cycle l2
          ty=atoms_cart(j)%type
          if(get_c_s(atoms(ty)%c_s) == "C") then
             if(k==0) then
                atty(1)=ty
                k=k+1
                m=k
             else
               l5: do m=1,k
                  if(ty==atty(m)) goto 1
               end do l5
               atty(k+1)=ty
               k=k+1
               m=k
             end if
1            at=m*one
             l3: do l=1,lenght
                if(j==cs_pair(1,l)) then
                   i1=i1+1
                   ishell(i1) = cs_pair(2,l)
                   st(i1)=at+0.01_r8_kind
                   exit l3
                elseif(j==cs_pair(2,l)) then
                   i1=i1+1
                   ishell(i1) = cs_pair(1,l)
                   st(i1)=at+0.01_r8_kind
                   exit l3
                end if
             end do l3
             write(file_gxcv,'(f5.2,3f15.7)') at,atoms_cart(i)%r(1:3)
          else if(get_c_s(atoms(j)%c_s) == "S") then
             cycle l1
          end if
       end do l2
    end do l1

    do i=1,lenght
       j = ishell(i)
       at=st(i)
       l4: do l=1,lenght
          if(j==cs_pair(1,l)) then
             k=cs_pair(2,l)
             exit l4
          elseif(j==cs_pair(2,l)) then
             k=cs_pair(1,l)
             exit l4
          end if
       end do l4
       x=atoms_frac(j)%r(1)
       y=atoms_frac(j)%r(2)
       z=atoms_frac(j)%r(3)
       if(abs(atoms_frac(j)%r(1)-atoms_frac(k)%r(1)) >= half) then
          if(atoms_frac(k)%r(1) > half) x=x+one
          if(atoms_frac(k)%r(1) <= half) x=x-one
       end if
       if(abs(atoms_frac(j)%r(2)-atoms_frac(k)%r(2)) >= half) then
          if(atoms_frac(k)%r(2) > half) y=y+one
          if(atoms_frac(k)%r(2) <= half) y=y-one
       end if
       if(lattice_calc) then
          if(abs(atoms_frac(j)%r(3)-atoms_frac(k)%r(3)) >= half) then
             if(atoms_frac(k)%r(3) > half) z=z+one
             if(atoms_frac(k)%r(3) <= half) z=z-one
          end if
       end if
       if(slab_calc) then
          x1=vect_s%v1(1)*x+vect_s%v2(1)*y
          y1=vect_s%v1(2)*x+vect_s%v2(2)*y
          z1=z
       elseif(lattice_calc) then
          x1=vect%v1(1)*x+vect%v2(1)*y+vect%v3(1)*z
          y1=vect%v1(2)*x+vect%v2(2)*y+vect%v3(2)*z
          z1=vect%v1(3)*x+vect%v2(3)*y+vect%v3(3)*z
      end if
       write(file_gxcv,'(f5.2,3f15.7)') at,x1,y1,z1
    end do

    call close_file_device(file_gxcv)

    if (allocated(ishell)) then
       deallocate (ishell, st, stat=status)
       if(status /= 0) call error_handler("MolMech: failed ISHELL deallocation")
    end if


  end subroutine epe_interface
  !******************************************************************

  !******************************************************************
  function name2type(at_nm,c_s)

    integer(kind=i4_kind) :: name2type
    character(len=len_name) :: at_nm
    character(len=1) :: c_s

    integer(kind=i4_kind) :: i

    name2type=0_i4_kind
    do i=1,n_species_types
       if(trim(at_nm) == trim(atoms(i)%name) .and. &
            c_s==get_c_s(atoms(i)%c_s)) then
          name2type=i
          exit
       end if
    end do

  end function name2type
  !******************************************************************

  !******************************************************************
  function get_c_s(num)
    
    character(len=1) :: get_c_s
    integer(kind=i4_kind) :: num
    
    if(num==0) get_c_s="C"
    if(num==1) get_c_s="S"

  end function get_c_s
  !******************************************************************

  !******************************************************************
  subroutine species_resorting()

    real(kind=r8_kind) :: rmax(3),rmin(3),dR(3)
    type(atom_coord) :: current_coor
    integer(kind=i4_kind) :: i,j,k,l,m,n,n1

    !search maximum and minimum values along X, Y and Z axes
    do i=1,n_species
       if(atoms_cart(i)%r(1) < rmin(1)) rmin(1)=atoms_cart(i)%r(1)
       if(atoms_cart(i)%r(2) < rmin(2)) rmin(2)=atoms_cart(i)%r(2)
       if(atoms_cart(i)%r(3) < rmin(3)) rmin(3)=atoms_cart(i)%r(3)
       if(atoms_cart(i)%r(1) > rmax(1)) rmax(1)=atoms_cart(i)%r(1)
       if(atoms_cart(i)%r(2) > rmax(2)) rmax(2)=atoms_cart(i)%r(2)
       if(atoms_cart(i)%r(3) > rmax(3)) rmax(3)=atoms_cart(i)%r(3)
    end do

    dR=abs(rmax-rmin)
    resort_axis=maxloc(dR)
    ! if resort_axis=1 resorting atoms along X axis
    ! if resort_axis=2 resorting atoms along Y axis
    ! if resort_axis=3 resorting atoms along Z axis
 
    ! resorting atoms (species)
    k=resort_axis(1)
    first_cycle: do i=2,n_species
       second_cycle: do j=1,i-1
          current_coor=atoms_cart(i)
          if(current_coor%r(k) > atoms_cart(j)%r(k)) cycle second_cycle
          if(current_coor%r(k) == atoms_cart(j)%r(k)) then
             do n=1,2
                n1=k+n
                if(n1 > 3) n1=n1-3
                if(current_coor%r(n1) > atoms_cart(j)%r(n1)) cycle second_cycle
                if(current_coor%r(n1) < atoms_cart(j)%r(n1)) exit
             end do
          end if
          third_cycle: do l=i-1,j,-1
             atoms_cart(l+1)=atoms_cart(l)
          end do third_cycle
          atoms_cart(j)=current_coor
          exit second_cycle
       end do second_cycle
    end do first_cycle

    do i=1,n_species
       n=atoms_cart(i)%neigh_sp
       do j=1,n_species
          m=atoms_cart(j)%initial_number
          if(n==m) then
             atoms_cart(i)%neigh_sp=j
             exit
          end if
       end do
    end do

print*,'============'
do i=1,n_species
print*,i,atoms_cart(i)%initial_number ,atoms_cart(i)%neigh_sp
end do
print*,'============'

  end subroutine species_resorting
  !******************************************************************

  !******************************************************************
  subroutine shutdown_species_data()

    integer(i4_kind) :: status,i

    deallocate(atoms_cart,stat=status)
    if(status /= 0) call error_handler("MolMech: failed ATOMS_CART deallocation")
    if(allocated(atoms_frac)) then
       deallocate(atoms_frac,stat=status)
       if(status /= 0) call error_handler("MolMech: failed ATOMS_FRAC deallocation")
    endif
    deallocate(atoms,stat=status)
    if(status /= 0) call error_handler("MolMech: failed ATOMS deallocation")
    if(allocated(im_coor)) then
       do i=1,n_species
          deallocate(im_coor(i)%r,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed IM_COOR%R deallocation(1)")
       end do
       deallocate(im_coor,stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR deallocation")
    end if
    if(allocated(im_coor_solv)) then
       do i=1,n_species
          deallocate(im_coor_solv(i)%r,stat=status)
          if(status /= 0) call error_handler( &
               "MolMech: failed IM_COOR%R deallocation(1)")
       end do
       deallocate(im_coor_solv,stat=status)
       if(status /= 0) call error_handler( &
            "MolMech: failed IM_COOR deallocation")
    end if

  end subroutine shutdown_species_data
  !******************************************************************

  !******************************************************************
  subroutine send_receive_species()

    integer(i4_kind) :: info,status,i

    if(comm_i_am_master()) then
       call comm_init_send(comm_all_other_hosts,msgtag_mm_send_species)

       call commpack(lattice_calc,info)
       if( info /= 0) call error_handler &
            ("send_receive_species: lattice_calc pack failed")
       call commpack(slab_calc,info)
       if( info /= 0) call error_handler &
            ("send_receive_species: slab_calc pack failed")
       call commpack(n_species,info)
       if( info /= 0) call error_handler &
            ("send_receive_species: n_species pack failed")
       do i=1,n_species
          call commpack(atoms_cart(i)%type,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms_cart(i)%type pack failed")
          call commpack(atoms_cart(i)%r(1),3,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms_cart(i)%r pack failed")
       end do

       call commpack(n_species_types,info)
       do i=1,n_species_types
          call commpack(atoms(i)%c_s,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms(i)%c_s pack failed")
          call commpack(atoms(i)%charge,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms(i)%charge pack failed")
       end do
       call comm_send()
    else
       call communpack(lattice_calc,info)
       if( info /= 0) call error_handler &
            ("send_receive_species: lattice_calc unpack failed")
       call communpack(slab_calc,info)
       if( info /= 0) call error_handler &
            ("send_receive_species: slab_calc unpack failed")
       call communpack(n_species,info)
       if( info /= 0) call error_handler &
            ("send_receive_species: n_species unpack failed")
       allocate(atoms_cart(n_species),stat=status)
       if(status /= 0) call error_handler &
               ("send_receive_species: allocation atoms_cart failed")
       do i=1,n_species
          call communpack(atoms_cart(i)%type,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms_cart(i)%type unpack failed")
          call communpack(atoms_cart(i)%r(1),3,1,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms_cart(i)%r unpack failed")
       end do

       call communpack(n_species_types,info)
       allocate(atoms(n_species_types),stat=status)
       if(status /= 0) call error_handler &
            ("send_receive_species: allocation atoms failed")
       do i=1,n_species_types
          call communpack(atoms(i)%c_s,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms(i)%c_s unpack failed")
          call communpack(atoms(i)%charge,info)
          if( info /= 0) call error_handler &
               ("send_receive_species: atoms(i)%charge unpack failed")
       end do
    end if

  end subroutine send_receive_species
  !******************************************************************

  !******************************************************************
  subroutine shutdown_species_on_slaves()

    integer(i4_kind) :: status

    deallocate(atoms_cart,stat=status)
    if(status /= 0) call error_handler &
         ("shutdown_species_on_slaves: deallocation atoms_cart failed")
    deallocate(atoms,stat=status)
    if(status /= 0) call error_handler &
         ("shutdown_species_on_slaves: deallocation atoms failed")

  end subroutine shutdown_species_on_slaves
  !******************************************************************

  !******************************************************************
end module species_module










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
module ewald_solv_module
!== Interrupt end of public interface of module =====              

!----modules used ------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  use common_data_module
  use tasks_main_options_module
  use species_module
  use cavity_module
  use ewald_module, only: erf_erfc
  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private

!------------ public functions and subroutines ------------------
  public calc_number_of_3D_images, init_ewald_3D_solv, calc_gauss_solv
  public shutdown_ewald_3D_solv, ewald_3d_contrib, ewald_3d_contrib_grad
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  integer(i4_kind) :: na,nb,nc
  real(r8_kind) :: ewald_param_solv,rcut_ew_solv,kcut_ew_solv,sqrt_pi
  integer(i4_kind), parameter :: max_kv = 200
  integer(i4_kind), allocatable :: k_indexes(:)
  integer(i4_kind) :: n_kvec
  real(r8_kind), allocatable :: gauss_multiplier(:)
  real(r8_kind), allocatable :: kvector(:,:)

contains
  !********************************************************
  subroutine calc_number_of_3D_images()
    !calculate coordinates of the 3D images of the unique cavity points
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind) :: amax,bmax,cmax,a_sin,b_sin,g_sin
    !------------ Executable code --------------------------------

    a_sin=sin(pi*cel%alpha/pi_degree)
    b_sin=sin(pi*cel%beta/pi_degree)
    g_sin=sin(pi*cel%gamma/pi_degree)

    amax=cel%a*b_sin*g_sin*half
    bmax=cel%b*g_sin*a_sin*half
    cmax=cel%c*a_sin*b_sin*half

    na=ceiling(rcut_ew_solv/amax); nb=ceiling(rcut_ew_solv/bmax); nc=ceiling(rcut_ew_solv/cmax)
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

  end subroutine calc_number_of_3D_images
  !********************************************************

  !********************************************************
  subroutine init_ewald_3D_solv()
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind) :: f
    real(r8_kind), parameter :: w_factor=one
    real(r8_kind), parameter :: accuracy=small
    !------------ Executable code --------------------------------

    ewald_param_solv=sqrt(((n_size*w_factor*pi**3)/volume**2)**(one/three))
    sqrt_pi=sqrt(pi)
    f=sqrt(-log(accuracy))
    rcut_ew_solv=f/ewald_param_solv
    kcut_ew_solv=two*f*ewald_param_solv

    call calc_kindexes_3D_solv()

  end subroutine init_ewald_3D_solv
  !********************************************************

  !********************************************************
  subroutine calc_kindexes_3D_solv()
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: max1,max2,max3
    real(r8_kind) :: drk,rk1(3),rk2(3),rk3(3),rk(3),cut2
    integer(i4_kind) :: i,j,k,status

    type store
       integer(i4_kind) :: index_store
       type(store), pointer :: next_data
    end type store
    type(store), target :: first_data
    type(store), pointer ::   current_data, tmp_data, del
    !------------ Executable code --------------------------------

    rk1=kvect%v1
    drk=sqrt(dot_product(rk1,rk1))
    max1=int(kcut_ew_solv/drk)+1

    rk2=kvect%v2
    drk=sqrt(dot_product(rk2,rk2))
    max2=int(kcut_ew_solv/drk)+1

    rk3=kvect%v3
    drk=sqrt(dot_product(rk3,rk3))
    max3=int(kcut_ew_solv/drk)+1

    if(max1 > max_kv .or. max2 > max_kv .or. max3 > max_kv) &
         call error_handler("MolMech-Solvation: Number of K vectors is too large")

    current_data=>first_data
    nullify(current_data%next_data)

    cut2=kcut_ew_solv*kcut_ew_solv

    n_kvec=0
    do i=0,max1
       rk1=i*kvect%v1
       do j=-max2,max2
          rk2=j*kvect%v2
          do k=-max3,max3
             if(i==0.and.j==0.and.k==0) cycle
             rk3=k*kvect%v3
             rk=rk1+rk2+rk3
             drk=dot_product(rk,rk)
             if(drk > cut2) cycle
             n_kvec=n_kvec+1
             allocate(tmp_data,stat=status)
             ASSERT(status == 0)
             tmp_data%index_store=1000000*(i+max_kv)+1000*(j+max_kv)+(k+max_kv)
             nullify(tmp_data%next_data)
             current_data%next_data =>tmp_data
             current_data => tmp_data
          end do
       end do
    end do

    allocate(k_indexes(n_kvec),stat=status)
    ASSERT(status == 0)
    current_data=>first_data
    do i=1,n_kvec
       del=>current_data
       current_data => current_data%next_data
       k_indexes(i)=current_data%index_store
       if(i > 1) deallocate(del)
    end do
    nullify(current_data)

  end subroutine calc_kindexes_3D_solv
  !********************************************************

  !********************************************************
  subroutine calc_gauss_solv()
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: ii,i,j,k,k_ind,status
    real(r8_kind) :: kv(3),dkv,ep2,factor
    !------------ Executable code --------------------------------

    allocate(gauss_multiplier(n_kvec),kvector(3,n_kvec),stat=status)
    ASSERT(status == 0)

    ep2=ewald_param_solv*ewald_param_solv

    do ii=1,n_kvec
       factor=two
       k_ind=k_indexes(ii)
       i=k_ind/1000000-max_kv
       j=(k_ind-(i+max_kv)*1000000)/1000-max_kv
       k=k_ind-(i+max_kv)*1000000-(j+max_kv)*1000-max_kv
       if(i==0) factor=one

       kv=i*kvect%v1+j*kvect%v2+k*kvect%v3
       kvector(:,ii)=kv
       dkv=dot_product(kv,kv)

       gauss_multiplier(ii)=factor*four*pi*exp(-dkv/(four*ep2))/(dkv*volume)
    end do
  end subroutine calc_gauss_solv
  !********************************************************

  !******************************************************
  function ewald_3d_contrib(rij) result(Aij)
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: Aij
    real(r8_kind) :: rij(3)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind) :: drij, drij_2,er,er2,erfunc,r(3),krij
    integer(i4_kind) :: i,j,k
    !------------ Executable code --------------------------------

    Aij=zero

    if(minimal_image) call image(rij)
    drij_2=dot_product(rij,rij); drij=sqrt(drij_2)

    if(drij > rcut_ew_solv) goto 1

    !Direct space
    er=ewald_param_solv*drij; er2=er*er
    erfunc=erf_erfc(er,ierfc)

    Aij=Aij+erfunc/drij

    if(.not. minimal_image) then
       do i=-na,na
          do j=-nb,nb
             do k=-nc,nc
                if(na==0 .and. nb==0 .and. nc==0) cycle
                r=rij+i*vect%v1+j*vect%v2+k*vect%v3
                drij_2=dot_product(r,r); drij=sqrt(drij_2)
                if(drij > rcut_ew_solv) cycle

                er=ewald_param_solv*drij; er2=er*er
                erfunc=erf_erfc(er,ierfc)

                Aij=Aij+erfunc/drij
             end do
          end do
       end do
    end if

1   continue
    !Reciprocal space
    do k=1,n_kvec
       krij=dot_product(kvector(:,k),rij)

       Aij=Aij+cos(krij)*gauss_multiplier(k)
    end do

  end function ewald_3d_contrib
  !******************************************************

  !******************************************************
  function ewald_3d_contrib_grad(rij,d_rij) result(dAij)
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    real(r8_kind) :: dAij(3)
    real(r8_kind) :: rij(3),d_rij(3,3)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(r8_kind) :: drij, drij_2,er,er2,erfunc,r(3),krij
    integer(i4_kind) :: i,j,k
    !------------ Executable code --------------------------------

    dAij=zero

    if(minimal_image) call image(rij)
    drij_2=dot_product(rij,rij); drij=sqrt(drij_2)

    if(drij > rcut_ew_solv) goto 1

    !Direct space
    er=ewald_param_solv*drij; er2=er*er
    erfunc=erf_erfc(er,ierfc)


    dAij=dAij+ &
         (-erfunc/drij_2-two*ewald_param_solv*exp(-er2)/(drij*sqrt_pi))*(rij/drij)

    if(.not. minimal_image) then
       do i=-na,na
          do j=-nb,nb
             do k=-nc,nc
                if(na==0 .and. nb==0 .and. nc==0) cycle
                r=rij+i*vect%v1+j*vect%v2+k*vect%v3
                drij_2=dot_product(r,r); drij=sqrt(drij_2)
                if(drij > rcut_ew_solv) cycle

                er=ewald_param_solv*drij; er2=er*er
                erfunc=erf_erfc(er,ierfc)

                dAij=dAij+ &
                     (-erfunc/drij_2-two*ewald_param_solv*exp(-er2)/(drij*sqrt_pi))*(r/drij)
             end do
          end do
       end do
    end if

1   continue
    !Reciprocal space
    do k=1,n_kvec
       krij=dot_product(kvector(:,k),rij)

       dAij=dAij-sin(krij)*gauss_multiplier(k)*kvector(:,k)
    end do

    dAij=matmul(dAij,d_rij)

  end function ewald_3d_contrib_grad
  !******************************************************

  !********************************************************
  subroutine shutdown_ewald_3D_solv()
    !------------ Modules used -----------------------------------
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: status
    !------------ Executable code --------------------------------

    if(allocated(k_indexes)) then
       deallocate(k_indexes,stat=status)
       ASSERT(status == 0)
    end if

    if(allocated(gauss_multiplier)) then
       deallocate(gauss_multiplier,kvector,stat=status)
       ASSERT(status == 0)
    end if

  end subroutine shutdown_ewald_3D_solv
  !********************************************************

end module ewald_solv_module

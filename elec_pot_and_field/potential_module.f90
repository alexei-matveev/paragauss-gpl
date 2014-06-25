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
!===============================================================
! Public interface of module
!===============================================================
module potential_module
!
!  Calculate electrostatic potential (nuclear and electronic parts)
!
!  A resulted value is a sum of potentials in equivalent points.
!  To get real magnitude of potential a calculated potential has
!  to be devided in a number of equivalent points
!
!  Author: AS
!  Date: 11/99
!
!================================================================
!== Interrupt of public interface of module =========
!================================================================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   2/01
! Description: potential due to point charges has been added
!
! Modification (Please copy before editing)
! Author: MF
! Date:   7/2000
! Description: adaptions to VPP
!
!----------------------------------------------------------------
!== Interrupt end of public interface of module =====

# include "def.h"
  use type_module          ! contains standard data types
  use datatype    ! user defined types
  use options_module
  use msgtag_module
  use comm_module
  use symmetry_data_module  ! provide ssym
  implicit none

  private
  save

!------------ public functions and subroutines ------------------

  public :: send_recv_space_point  ! ()

  public start_read_poten_e,read_poten_e_3,get_poten_n,get_poten_pc, &
       dealloc_space_points, deallocate_pot, &
       bounds_calc_poten, get_bounds_poten, &
       bounds_free_poten,send_receive_poten, fill_points, &
       destroy_poten_file , &
       poten_integral_open,poten_integral_close

  type, public :: spacepoint_type
     integer (i4_kind) :: N_equal_points
     ! Number of partners of unique point charge
     real (r8_kind), allocatable :: position(:, :)
     ! positions of partners ("equal points") (3 x N_equal_points)
  end type spacepoint_type

  ! FIXME: %lower1 and %upper1 are obscure and may be removed,
  !        see bounds_module.f90. After that this struct is just
  !        a plain array.
  type,public :: poten_bounds
     integer(kind=i4_kind)                :: lower1,upper1
     ! start bounds
     integer(kind=i4_kind),pointer        :: item_arr(:)
  end type poten_bounds

  ! These two are not protected:
  type (spacepoint_type), allocatable, target, public :: point_in_space(:)
  integer (i4_kind), public :: N_points

  ! Potential of  electronic or nuclear  charge at the  unique surface
  ! points of the cavity, respectively.
  real (r8_kind), allocatable, public, protected :: V_pot_e(:)
  real (r8_kind), allocatable, public, protected :: V_pot_n(:)
  real (r8_kind), allocatable, public, protected :: V_pot_pc(:)

#ifdef WITH_EFP
  real (r8_kind), allocatable, public :: V_pot_mp(:)  ! not protected
  real (r8_kind), allocatable, public :: V_pot_id(:)  ! not protected
  real (r8_kind), allocatable, public :: V_pot_id1(:) ! not protected
#endif

  type(poten_bounds) :: bounds
  real(kind=r8_kind), allocatable :: V_pot_e_tmp(:)

  integer(kind=i4_kind) :: n_irrep

!================================================================
! End of public interface of module
!================================================================
!To calculate electrostatic potential:
!1. Define points where electrostatic potential has to be calculated P(3,N)
!2. call fill_points(P,N)
!3. call potential_calculate() - integrals
!4. call main_scf() or use saved density matrix (main_scf can be called
!   in place before calculation of electrostatic field)
!5. call start_read_poten_e()
!6. call get_poten_n()
!7. call get_poten_pc()
!8.Here you can use V_pot_e, V_pot_n, V_pot_pc (please remember that
!  the real potential is equal to V_pot_.../N_equal_points)
!9. call deallocate_pot()
!10.call dealloc_space_points()
!11.call bounds_free_poten()
!12.call post_dens_master()
!10. If you want to keep integrals in memory, do not forget to deallocate them in the end
!           if ( .not. options_integrals_on_file() ) then !!!!!!!!!!!!!!!!!
!              call integralstore_deallocate_pcm()
!              if ( comm_parallel() ) then
!                 call comm_init_send(comm_all_other_hosts,msgtag_intstore_dealloc)
!                 call comm_send()
!              endif
!            end if

contains

  !******************************************************
  subroutine fill_points (point_array, N_total)
    use group_module, only : ylm_trafos,sub_group,group_coset, &
         symm_transformation_int,group_num_el,group_coset_decomp
    implicit none
    real (r8_kind), intent (in) :: point_array(3, N_total)
    integer (i4_kind), intent (in) :: N_total
    !== End of interface ==========================================

    integer(kind=i4_kind) :: n_equal,n_equal_check,n_size
    real(kind=r8_kind)    :: position(3),position2(3),position3(3)
    type(sub_group) :: local_groups
    type(group_coset) :: cosets
    type(symm_transformation_int) :: point_trafos
    logical, allocatable :: help_dim(:)
    type(spacepoint_type),allocatable :: buffer(:)

    real(kind=r8_kind),parameter :: small = 1.0e-10_r8_kind
    integer(kind=i4_kind) :: i,n,j,k,l,status
    !------------ Executable code ---------------------------------

    allocate (help_dim(N_total), stat=status)
    if ( status /= 0) call error_handler( &
         "potential_module: allocation of help_dim is failed")
    help_dim=.false.

    allocate(buffer(N_total),stat=status)
    if ( status /= 0) call error_handler( &
         "potential_module: allocation of buffer is failed")

    n_size=0
    l=0
    i_N_tot: do i=1,N_total
       if(help_dim(i)) cycle i_N_tot
       n_size=n_size+1
       ! reorder coordinates of points as
       ! (x,y,z) --> (x,z,y) in order to comply with the
       ! convention for angular momentum l=1
       position(1) = point_array(1, i)
       position(2) = point_array(3, i)
       position(3) = point_array(2, i)
       !
       ! determine local symmetry groups
       !
       ! now apply all symmetry operations to the position of the
       ! surface point
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
          endif
       enddo


       ! allocate group elements
       local_groups%num_el = n
       allocate(local_groups%elements(n))
       if ( status /= 0) call error_handler( &
            "potential_module: allocation LOCAL_GROUPS  is failed")

       ! fill up group elements
       n = 0
       do j=1,group_num_el
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,j),position)
          if (dot_product(position2-position,position2-position) < small) then
             n = n+1
             local_groups%elements(n) = j
          end if
       enddo
       !
       ! now determine symmetry equivalent atoms
       !
       call group_coset_decomp(n_equal,local_groups,&
            cosets,point_trafos%matrix)
       !
       ! search of positions of equal atoms
       !
       allocate(buffer(n_size)%position(3,n_equal), stat=status)
       if (status .ne. 0) call error_handler( &
            "elec_static_field_module: allocation of BUFFER%POSITION  failed")

       n_equal_check=0
       j_n_equal: do j=1,n_equal
          position2 = MATMUL(ylm_trafos(1)%matrix(:,:,cosets%elements(1,j)),position)
          position3(1) = position2(1)
          position3(2) = position2(3)
          position3(3) = position2(2)
          k_N_tot: do k=1,N_total
             if(.not.help_dim(k)) then
                if (sqrt(dot_product(position3-point_array(:,k), &
                     position3-point_array(:,k))) <= small) then

                   help_dim(k)=.true.
                   n_equal_check=n_equal_check+1
                   buffer(n_size)%position(:,j)=point_array(:,k)
                   buffer(n_size)%n_equal_points=n_equal
                endif
             endif
          end do k_N_tot
       end do j_n_equal
       if (n_equal_check /= n_equal) call error_handler( &
            "potential_module: The point array to calculate electrostatic &
            & potential is not corresponded to the chosen symmetry")

       deallocate(local_groups%elements,point_trafos%matrix,cosets%elements,stat=status)
       if (status .ne. 0 ) call error_handler( &
            "potential_module: deallocation of points HELPERS is failed")
    end do i_N_tot

    N_points=n_size
    allocate(point_in_space(N_points),stat=status)
    if ( status /= 0) call error_handler( &
         "potential_module: allocation POINTS_IN_SPACE is  failed")
    do i=1,n_size
       point_in_space(i)%n_equal_points=buffer(i)%n_equal_points

       allocate(point_in_space(i)%position(3,buffer(i)%n_equal_points), stat=status)
       if (status .ne. 0 ) call error_handler( &
            "potential_module: allocation of POINTS_IN_SPACE(i)%position is failed")

       point_in_space(i)%position=buffer(i)%position

       deallocate(buffer(i)%position, stat=status)
       if (status .ne. 0 ) call error_handler( &
            "potential_module: deallocation of BUFFER(i)%position is failed")
    enddo
    deallocate(buffer,stat=status)
    if ( status /= 0) call error_handler( &
         "potential_module: deallocation of BUFFER is failed")
    deallocate(help_dim,stat=status)
    if ( status /= 0) call error_handler( &
         "potential_module: deallocation of HELP_DIM  is failed")

    if (comm_parallel()) call send_space_point()

  end subroutine fill_points
  !********************** ********************************

  subroutine send_recv_space_point ()
    use comm, only: comm_rank
    use comm_module, only: comm_save_recv, comm_master_host
    use msgtag_module, only: msgtag_space_point
    implicit none
    ! *** end of interface ***

    if (comm_rank () == 0) then
       call send_space_point ()
    else
       call comm_save_recv (comm_master_host, msgtag_space_point)
       call receive_space_point ()
    endif
  end subroutine send_recv_space_point

  !******************************************************
  subroutine send_space_point ()
    ! send representative surface points of the cavity
    !** End of interface *****************************************

    integer (i4_kind) :: i, length, status
    type (spacepoint_type), pointer :: ps

    call comm_init_send (comm_all_other_hosts, msgtag_space_point)

    call commpack(N_points,status)
    do i=1,N_points
       ps => point_in_space(i)
       call commpack(ps%N_equal_points,status)
       length = ps%N_equal_points*3
       call commpack(ps%position(1,1), length, 1, status)
    enddo
    call comm_send()
  end subroutine send_space_point
  !******************************************************

  !******************************************************
  subroutine receive_space_point ()
    !
    ! Receive representative surface points of the cavity
    !
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind)           :: i, length, status

    call communpack(N_points,status)
    allocate(point_in_space(N_points),stat=status)
    if ( status .ne. 0) call error_handler( &
         "potential_module:receive_space_point: allocated point_in_space failed" )

    do i=1,N_points
       call communpack(point_in_space(i)%N_equal_points,status)
       allocate(point_in_space(i)%position(3,point_in_space(i)%N_equal_points),stat=status)
       if ( status .ne. 0) call error_handler( &
            "potential_module:receive_space_point: allocated point_in_space%position failed" )
       length = point_in_space(i)%N_equal_points*3
       call communpack(point_in_space(i)%position(1,1), length, 1, status)
    enddo
  end subroutine receive_space_point
  !******************************************************

  !******************************************************
  subroutine dealloc_space_points
    ! deallocate "point_in_space" (information about surface points
    ! of cavity
    !** End of interface *****************************************

    integer(kind=i4_kind)           :: i, status

    do i=1,N_points
       deallocate(point_in_space(i)%position,stat=status)
       if ( status .ne. 0) call error_handler( &
            "potential_module: deallocated point_in_space%position failed" )
    enddo

    ! With only allocatable components this should be enough:
    deallocate (point_in_space, stat=status)
    if ( status .ne. 0) call error_handler( &
         "potential_module: deallocated point_in_space failed" )
  end subroutine dealloc_space_points
  !******************************************************

  !*********************************************************
  subroutine start_read_poten_e
   !** End of interface *****************************************

    if(comm_parallel()) then
       call comm_init_send(comm_all_other_hosts,msgtag_start_poten)
       call comm_send()
    endif

    call read_poten_e_3
    call send_receive_poten

  end subroutine start_read_poten_e
  !*********************************************************

  !*********************************************************
  subroutine read_poten_e_3
    use density_data_module,only: densmat
    use readwriteblocked_module
    use integralstore_module, only: integralstore_3c_poten
    implicit none
    !** End of interface *****************************************

    type(readwriteblocked_tapehandle) :: th_poten
    real(kind=r8_kind),pointer :: poten_int(:)
    integer (i4_kind) :: item_arr_poten, my_ind, n_irrep, i_gamma, mul
    logical ::  spin_polarized,integrals_on_file
    integer (i4_kind) :: n, status, m, k, i_meta, i_last
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
#ifdef _VECTOR
    integer (i4_kind) :: i_mn, n_mn
    integer(kind=i4_kind), allocatable :: m_i(:),n_i(:)
#endif

#ifdef _VECTOR
    DPRINT '_e_3: compiled for vector machines'
#else
    DPRINT '_e_3: compiled for scalar machines'
#endif
    integrals_on_file=options_integrals_on_file()
    DPRINT 'potential_module:read_poten_e_3: on_file=',integrals_on_file
    spin_polarized = symmetry_data_n_spin() > 1

    n_irrep = symmetry_data_n_irreps()

    allocate(dim_irrep(n_irrep),stat=status)
    if ( status .ne. 0) call error_handler( &
         "potential_module:read_poten_e_3: allocated dim_irrep failed" )

    do n=1,n_irrep
       dim_irrep(n)  = symmetry_data_dimension(n)
    enddo

    my_ind=comm_myindex()

    item_arr_poten=bounds%item_arr(my_ind)

    if (item_arr_poten == 0) then
       deallocate(dim_irrep,stat=status)
       if ( status .ne. 0) call error_handler( &
            "potential_module:read_poten_e_3: deallocated dim_irrep failed" )
       return
    endif

    allocate(V_pot_e_tmp(item_arr_poten), stat=status)
    if ( status .ne. 0) call error_handler( &
         "potential_module:read_poten_e_3: allocated  V_pot_e_tmp failed" )
    V_pot_e_tmp=0.0_r8_kind

    if ( integrals_on_file ) then
       ! open the integral files ----------------------
       call poten_integral_open(th_poten)
    else
         i_meta=1
    endif

#ifndef _VECTOR
    if (integrals_on_file) then
       allocate( poten_int(item_arr_poten),STAT=status)
       if(status.ne.0) call error_handler&
            ("potential_module: allocation of poten_int failed")
    endif
#endif

    do i_gamma = 1,n_irrep
#ifndef _VECTOR
      do m=1,dim_irrep(i_gamma)
         do n=1,m
            mul = 2
            if (n == m) mul = 1
            ! read in integrals file
            if (integrals_on_file) then
               call readwriteblocked_read(poten_int,th_poten)
            else
               i_last=i_meta+item_arr_poten-1
               poten_int=> integralstore_3c_poten(i_meta:i_last)
               i_meta=i_last+1
            endif
            do k=1,item_arr_poten
               if(spin_polarized) then
                  V_pot_e_tmp(k)=V_pot_e_tmp(k)-poten_int(k)* &
                     mul * (densmat(i_gamma) % m(n, m, 1) + densmat(i_gamma) % m(n, m, 2))
               else
                  V_pot_e_tmp(k)=V_pot_e_tmp(k)-poten_int(k)* &
                      mul * densmat(i_gamma) % m(n, m, 1)
               endif
            enddo
         enddo
      enddo

#else

      n_mn = ( dim_irrep(i_gamma) * (dim_irrep(i_gamma) + 1) ) / 2
      allocate(m_i(n_mn),n_i(n_mn))
      i_mn=1
      do m=1,dim_irrep(i_gamma)
         do n=1,m
           m_i(i_mn)=m
           n_i(i_mn)=n
           i_mn=i_mn+1
         enddo
      enddo
      if ( integrals_on_file ) then
            allocate( poten_int(item_arr_poten*n_mn),STAT=status )
            if(status.ne.0) call error_handler &
                 ("potential_module: allocation of poten_int failed")
            call readwriteblocked_read(poten_int,th_poten)
               do k=1,item_arr_poten
                 do i_mn = 1, n_mn
                   mul = 2
                   if (n_i(i_mn) == m_i(i_mn)) mul = 1
                   if(spin_polarized) then
                    V_pot_e_tmp(k)=V_pot_e_tmp(k)-&
                       poten_int(k+(i_mn-1)*item_arr_poten)* &
                       mul * (densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1)&
                       + densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 2))
                   else
                    V_pot_e_tmp(k)=V_pot_e_tmp(k)-&
                        poten_int(k+(i_mn-1)*item_arr_poten)* &
                        mul * densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1)
                   endif
                 enddo !i_mn
               enddo !k
            deallocate(poten_int,STAT=status)
            if(status.ne.0) call error_handler &
                 ("potential_modul: deallocation of poten_int failed")
      else ! integrals not on file
               do k=1,item_arr_poten
                 do i_mn = 1, n_mn
                   mul = 2
                   if (n_i(i_mn) .eq. m_i(i_mn)) mul = 1
                   if(spin_polarized) then
                    V_pot_e_tmp(k)=V_pot_e_tmp(k)-&
                       integralstore_3c_poten(i_meta-1+k+(i_mn-1)*item_arr_poten)* &
                       mul * (densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1)&
                       + densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 2))
                   else
                    V_pot_e_tmp(k)=V_pot_e_tmp(k)-&
                        integralstore_3c_poten(i_meta-1+k+(i_mn-1)*item_arr_poten)* &
                        mul * densmat(i_gamma) % m(n_i(i_mn), m_i(i_mn), 1)
                   endif
                 enddo !i_mn
               enddo !k
               i_meta=i_meta+item_arr_poten*n_mn
      endif
      deallocate(m_i,n_i)
#endif

    enddo

#ifndef _VECTOR
      if (integrals_on_file) then
         deallocate(poten_int,STAT=status)
         if(status.ne.0) call error_handler&
              ("potential_modul: deallocation of poten_int failed")
      endif
#endif

    if ( integrals_on_file ) then
       ! close the integral files ----------------------
            call poten_integral_close(th_poten)
    endif

    deallocate(dim_irrep,stat=status)
    if ( status .ne. 0) call error_handler( &
         "potential_module:read_poten_e_3: deallocated dim_irrep failed(1)" )

  end subroutine read_poten_e_3
  !*********************************************************

  !*********************************************************
  subroutine send_receive_poten
    !** End of interface *****************************************

    integer(kind=i4_kind) :: my_ind,poten_item_arr,bound
    integer(kind=i4_kind) :: status,info,i_pr
    integer (i4_kind) :: i

    if(comm_parallel()) then
       if(comm_i_am_master()) then
          if(.not.allocated(V_pot_e)) then
             allocate(V_pot_e(N_points),stat=status)
             if ( status .ne. 0) call error_handler( &
                  "potential: allocated of V_pot_e failed" )
          endif

          i_pr=comm_get_n_processors()
          bound=0
          do i=1,i_pr
             if (i == 1) then
                poten_item_arr=bounds%item_arr(1)
                if(poten_item_arr == 0) cycle
                V_pot_e(bound+1:bound+poten_item_arr)= V_pot_e_tmp
                deallocate(V_pot_e_tmp,stat=status)
                if ( status .ne. 0) call error_handler( &
                     "potential: deallocated of V_pot_e_tmp failed(1)" )
             else
                call comm_save_recv(i,msgtag_finish_poten)

                call communpack(poten_item_arr,info)
                if( info.ne.0) call error_handler &
                     ("potential:unpack of poten_item_arr failed")
                if (poten_item_arr == 0) cycle

                call communpack(V_pot_e(bound+1:bound+poten_item_arr),poten_item_arr,1,info)
                if( info.ne.0) call error_handler &
                     ("potential:unpack of V_pot_e failed")
             endif
             bound=bound+poten_item_arr
          enddo
       else
          call comm_init_send(comm_master_host,msgtag_finish_poten)

          my_ind=comm_myindex()

          call commpack(bounds%item_arr(my_ind),info)
          if( info.ne.0) call error_handler &
               ("potential:pack of item_arr failed")

          if(bounds%item_arr(my_ind) /= 0) then
             call commpack(V_pot_e_tmp,bounds%item_arr(my_ind),1,info)
             if( info.ne.0) call error_handler &
                  ("potential:pack of V_pot_e_tmp failed")

             deallocate(V_pot_e_tmp,stat=status)
             if ( status .ne. 0) call error_handler( &
                  "potential: deallocated of V_pot_e_tmp failed(2)" )
          endif
          call comm_send()
       endif
    else
       if (.not.allocated(V_pot_e)) then
          allocate(V_pot_e(N_points),stat=status)
          if ( status .ne. 0) call error_handler( &
               "potential: allocated of V_pot_e failed(1)" )
       endif
       V_pot_e=V_pot_e_tmp

       deallocate(V_pot_e_tmp,stat=status)
       if ( status .ne. 0) call error_handler( &
            "potential: deallocated of V_pot_e_tmp failed(3)" )

    endif

  end subroutine send_receive_poten
  !*********************************************************

  !*********************************************************
  subroutine get_poten_n ()
    !
    ! Fills  public  protected  module  variable V_pot_n(:)  with  the
    ! electrostatic potential  of atomic cores. Note  that the nuclear
    ! charge of the ECPs is shielded by core electrons.
    !
    use unique_atom_module, only: unique_atoms, N_unique_atoms
    use unique_atom_methods, only: core_charge => unique_atom_core_charge
    implicit none
    !** End of interface *****************************************

    real (r8_kind) :: z, dist
    integer (i4_kind) :: na, eq_a, eq_b, n_equal_points, i, status
    real (r8_kind), pointer :: xa(:,:), xb(:,:)

    allocate (V_pot_n(N_points), stat=status)
    if ( status .ne. 0) call error_handler( &
         "potential:get_poten_n: allocated V_pot_n failed" )

    V_pot_n = 0.0
    do i = 1, N_points
       xb => point_in_space(i) % position
       n_equal_points = point_in_space(i) % N_equal_points
       unique_a: do na = 1, N_unique_atoms
          z = core_charge (unique_atoms(na)) ! z - zc for ECPs
          xa => unique_atoms(na) % position

          equal_a: do eq_a = 1, unique_atoms(na) % N_equal_atoms
             equal_p: do eq_b = 1, n_equal_points
                dist = sqrt (sum ((xa(:, eq_a) - xb(:, eq_b))**2))
                if (dist < 0.001_r8_kind) dist = 0.001_r8_kind
                V_pot_n(i) = V_pot_n(i) + z / dist
             enddo equal_p
          enddo equal_a
       enddo unique_a
    enddo
  end subroutine get_poten_n
  !*********************************************************

  !*********************************************************
  subroutine get_poten_pc ()
    use pointcharge_module, only : pointcharge_array, pointcharge_N
#ifdef WITH_EPE
    use ewaldpc_module, only : ewpc_array, EWPC_N
#endif
    implicit none
    ! *** end of interface ***

    real(kind=r8_kind)          :: z1,dist
    integer (i4_kind) :: na, eq_a, eq_b, n_equal_points, i, status
    real(kind=r8_kind),pointer  :: xa(:,:),xb(:,:)

    allocate(V_pot_pc(N_points),stat=status)
    if ( status .ne. 0) call error_handler( &
         "potential:get_poten_pc: allocated V_pot_pc failed" )
    V_pot_pc=0.0_r8_kind
    do i=1,N_points
       xb => point_in_space(i)%position
       n_equal_points = point_in_space(i)%N_equal_points
       unique_a: do na=1,pointcharge_N
          z1 = pointcharge_array(na)%Z
          xa => pointcharge_array(na)%position

          equal_a: do eq_a=1,pointcharge_array(na)%N_equal_charges
             equal_p: do eq_b=1,n_equal_points
                dist = sqrt(sum((xa(:,eq_a)-xb(:,eq_b))**2))
                if(dist < 0.001_r8_kind) dist=0.001_r8_kind
                V_pot_pc(i) = V_pot_pc(i) + z1/dist
             enddo equal_p
          enddo equal_a
       enddo unique_a

#ifdef WITH_EPE
       unique_ae: do na=1,ewpc_N
          z1 = ewpc_array(na)%Z
          xa => ewpc_array(na)%position

          equal_ae: do eq_a=1,ewpc_array(na)%N_equal_charges
             equal_pe: do eq_b=1,n_equal_points
                dist = sqrt(sum((xa(:,eq_a)-xb(:,eq_b))**2))
                if(dist < 0.001_r8_kind) dist=0.001_r8_kind
                V_pot_pc(i) = V_pot_pc(i) + z1/dist
             enddo equal_pe
          enddo equal_ae
       enddo unique_ae
#endif
    enddo
  end subroutine get_poten_pc
  !*********************************************************

  !*********************************************************
  subroutine deallocate_pot
    !** End of interface *****************************************

    integer(kind=i4_kind)       :: status

    if(allocated(V_pot_n)) then
       deallocate(V_pot_n,stat=status)
       if ( status .ne. 0) call error_handler( &
            "potential:deallocation of V_pot_n failed" )
    endif
    if(allocated(V_pot_e)) then
       deallocate(V_pot_e,stat=status)
       if ( status .ne. 0) call error_handler( &
            "potential:deallocation of V_pot_e failed" )
    endif
    if(allocated(V_pot_pc)) then
       deallocate(V_pot_pc,stat=status)
       if ( status .ne. 0) call error_handler( &
            "potential:deallocation of V_pot_pc failed" )
    endif
#ifdef WITH_EFP
    if(allocated(V_pot_mp)) then
       deallocate(V_pot_mp,stat=status)
       ASSERT(status==0)
    endif
#endif

  end subroutine deallocate_pot
  !*********************************************************

  !*********************************************************
  subroutine bounds_calc_poten ()
    !
    ! Fills global module var bounds.
    !
    use comm, only: comm_size, comm_same
    implicit none
    !** End of interface *****************************************

    real(kind=r8_kind)      :: rhelp
    integer(kind=i4_kind)   :: i,alloc_stat, i_pr

    ! If the input is the same, the output will be the same:
    ASSERT(comm_same(N_points))

    i_pr = comm_size ()

    allocate (bounds%item_arr(i_pr),STAT=alloc_stat)
    if(alloc_stat.ne.0) call error_handler&
         ("potential: allocation of BOUNDS failed ")
    bounds%item_arr=0

    if (comm_parallel()) then
       rhelp = real(N_points,kind=r8_kind)

          do i=1,i_pr
             bounds%item_arr(i) = ceiling(rhelp/i_pr)
             rhelp = rhelp - real(bounds%item_arr(i),kind=r8_kind)
             i_pr = i_pr - 1
          enddo

          bounds%lower1 = bounds%item_arr(1)+1
          bounds%upper1 = bounds%item_arr(1) + bounds%item_arr(2)

    else  ! no parallel processing

       bounds%item_arr(1) = N_points
       bounds%lower1 = 1
       bounds%upper1 = 0
    endif! end of if( comm_parallel )
  end subroutine bounds_calc_poten
  !******************************************************************

  !*************************************************************
  subroutine bounds_free_poten()
    !  Purpose: Release the private bounds pointer
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: status
    !------------ Executable code --------------------------------

    if (comm_parallel() .and. comm_i_am_master()) then
       call comm_init_send(comm_all_other_hosts,msgtag_free_bnds_ptn)
       call comm_send()
    endif

    deallocate (bounds%item_arr,STAT=status)
    if (status /= 0) call error_handler&
         ("potential: deallocation of BOUNDS failed ")
  end subroutine bounds_free_poten
  !*************************************************************

  !*************************************************************
  subroutine get_bounds_poten(bounds1)
    !** End of interface *****************************************
    !------------ Declaration of formal parameters ---------------
    type(poten_bounds),intent(out)          :: bounds1
    integer(kind=i4_kind)             :: alloc_stat,i_pr
    !------------ Executable code --------------------------------

    bounds1%lower1 = bounds%lower1
    bounds1%upper1 = bounds%upper1
    i_pr=comm_get_n_processors()
    allocate(bounds1%item_arr(i_pr),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("potential : allocation of BOUNDS1 failed ")
    bounds1%item_arr = bounds%item_arr

  end subroutine get_bounds_poten
  !***************************************************************

  !**********************************************************
  subroutine poten_integral_open(th_poten2)
    !  Purpose: decide which files to open, perform 'startread'
    !** End of interface *****************************************
    ! Modules:
    use filename_module, only: tmpfile
    use readwriteblocked_module
    !------------ Declaration of formal parameters ---------------
    type(readwriteblocked_tapehandle),intent(out) :: th_poten2
    !------------ Declaration of local variables -----------------
!!$ integer(kind=i4_kind)        :: alloc_stat
!!$ integer(kind=i4_kind)        :: mynumber
!!$ type(poten_bounds)             :: bounds2

    !------------ Executable code --------------------------------

!!$    mynumber = comm_myindex()
!!$
!!$    allocate(bounds2%item_arr(comm_get_n_processors()),STAT=alloc_stat)
!!$       if(alloc_stat.ne.0) then
!!$          call error_handler &
!!$               ("potential: allocation of BOUNDS2 failed")
!!$       endif
!!$
!!$    call get_bounds_poten(bounds2)
!!$
!!$    if( bounds2%item_arr(mynumber).eq.0) then
!!!!!!!       call error_handler( "prepare_integral_open: both bounds 0" )
!!$    else
       call readwriteblocked_startread(trim(tmpfile('poten.dat')), th_poten2)
!!$    endif

  end subroutine poten_integral_open
  !*************************************************************

  !*************************************************************
  subroutine poten_integral_close(th_poten1)
    !  Purpose: performs a 'stopread' on the files contained
    !           in the specified tapehandle,
    !           returns the corresponding io_unit and
    !           closes the files
    !** End of interface *****************************************
    use readwriteblocked_module
    !------------ Declaration of formal parameters ---------------
    type(readwriteblocked_tapehandle),intent(inout) :: th_poten1
    call readwriteblocked_stopread(th_poten1)
  end subroutine poten_integral_close
  !*************************************************************

  !*************************************************************
  subroutine destroy_poten_file
    ! Purpose: to delete files of potential matrix elements
    !          from disk
    !** End of interface *****************************************
    use filename_module, only: tmpfile
    use readwriteblocked_module

    logical :: yes
    type(readwriteblocked_tapehandle) :: th_poten

    inquire(file=trim(tmpfile('poten.dat')), exist=yes)
    if (yes) then
       call readwriteblocked_startread(trim(tmpfile('poten.dat')), th_poten)
       call readwriteblocked_stopread(th_poten,status='delete')
    endif

  end subroutine destroy_poten_file
  !*************************************************************

end module potential_module


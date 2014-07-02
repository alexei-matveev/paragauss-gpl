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
  !===================================================================
! Public interface of module
  !===================================================================
module print_module
!-------------- Module specification ---------------------------
!
!  This file contains all kinds of output routines.
!  Should be self-explanatory
!  TG,FN
!
! defines generic subroutine:
!
!     subroutine printout(unit,item,name)
!       use type_module
!       integer(kind=i4_kind),   intent(in)  :: unit
!       xxx,                     intent(in)  :: item
!       character*(*), optional, intent(in)  :: name
!     where type xxx of item can be:
!       integer(kind=i4_kind),dimension(:)
!       real(kind=r8_kind),pointer,dimension(:)
!       real(kind=r8_kind),pointer,dimension(:,:)
!       type(arrmat2),dimension(:)
!       type(arrmat3),dimension(:)
!       
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: TB
! Date:   1/96
! Description: included former file outroutines.f90
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

  implicit none
  save

  !== Interrupt end of public interface of module ====================
 
  interface printout
     module procedure print_real_vec
     module procedure print_real_mat
     module procedure print_int_vec
     module procedure print_arrmat3_vec
     module procedure print_arrmat2_vec
  end interface
        
!------------ public functions and subroutines ------------------
public banner, print_int_vec, print_real_vec, print_real_mat, &
       print_arrmat3_vec, print_arrmat2_vec


  !===================================================================
  ! End of public interface of module
  !===================================================================

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

subroutine banner(unit,s)
  ! prints s in banner form to unit
  use type_module
  implicit none 
  integer(kind=i4_kind), intent(in) :: unit
  character*(*), intent(in)         :: s
  !** End of interface *****************************************
  write(unit,*) '****************************************'
  write(unit,*)
  write(unit,*) s
  write(unit,*)
  write(unit,*) '****************************************'
end subroutine banner

subroutine print_int_vec(unit,item,name)
  use type_module
  implicit none
  integer(kind=i4_kind), intent(in)               :: unit
  integer(kind=i4_kind), dimension(:), intent(in) :: item
  character*(*),optional, intent(in)              :: name
  !** End of interface *****************************************
  integer(kind=i4_kind)                     :: i
  integer(kind=i4_kind)                     :: dimlo,dimup
  dimlo=lbound(item,1)
  dimup=ubound(item,1)
  write(unit,*)
  if (present(name)) then
     write(unit,*) 'Vektor    :',name
  endif
  write(unit,*) 'Dimension :',dimlo,' - ',dimup
  do i=dimlo,dimup
     write(unit,*) item(i)
  enddo
  write(unit,*)
end subroutine print_int_vec

subroutine print_real_vec(unit,item,name)
  use type_module
  implicit none
  integer(kind=i4_kind), intent(in)              :: unit
  real(kind=r8_kind), dimension(:), intent(in)   :: item
  character*(*),optional, intent(in)             :: name
  !** End of interface *****************************************
  integer(kind=i4_kind)                     :: i
  integer(kind=i4_kind)                     :: dimlo,dimup
  dimlo=lbound(item,1)
  dimup=ubound(item,1)
  write(unit,*)' '
   if (present(name)) then
     write(unit,*) 'Vektor    :',name
  endif

  write(unit,*) 'Dimension :',dimlo,' - ',dimup
  write(unit,*) ( item(i), i=dimlo,dimup)
  write(unit,*) ' '
end subroutine print_real_vec

subroutine print_real_mat(unit,item,name)
  use type_module
  implicit none
  integer(kind=i4_kind), intent(in)                    :: unit
  real(kind=r8_kind),dimension(:,:), intent(in)        :: item
  character*(*),optional, intent(in)                   :: name
  !** End of interface *****************************************
  integer(kind=i4_kind)  :: dimlo1,dimup1,dimlo2,dimup2,i,j
  dimlo1=lbound(item,1)
  dimup1=ubound(item,1)
  dimlo2=lbound(item,2)
  dimup2=ubound(item,2)
  write(unit,*)' '
  if (present(name)) then
     write(unit,*) 'Vektor    :',name
  endif

  write(unit,*) 'Dimension :',dimlo1,' - ',dimup1
  write(unit,*) '           ',dimlo2,' - ',dimup2
  do i=dimlo1,dimup1
     write(unit,*) ( item(i,j), j=dimlo2,dimup2)
  enddo
  write(unit,*) ' '
end subroutine print_real_mat

subroutine print_arrmat3_vec(unit,item,name)
  use type_module
  use datatype
  implicit none
  integer(kind=i4_kind), intent(in)     :: unit
  character*(*),optional, intent(in)    :: name
  type(arrmat3), intent(in)             :: item(:)
  !** End of interface *****************************************
  integer(kind=i4_kind)                 :: dimlo1,dimup1,i,j,k,l
  integer(kind=i4_kind),allocatable     :: dimi(:,:)
  dimlo1 = lbound(item,1)
  dimup1 = ubound(item,1)
  allocate( dimi(dimup1-dimlo1+1,3))
  do i=dimlo1,dimup1
     dimi(i,1) = ubound(item(i)%m,1)
     dimi(i,2) = ubound(item(i)%m,2)
     dimi(i,3) = ubound(item(i)%m,3)
  enddo
  write(unit,*)' '
  if (present(name)) then
     write(unit,*) 'Vektor    :',name
     write(unit,*)' Dimensionen:'
     write(unit,*)'IRREPS      :',dimlo1,' - ',dimup1
     do i=dimlo1,dimup1
        write(unit,*) dimi(i,1),' ',dimi(i,2),' ',dimi(i,3)
     enddo
  endif
  do i = dimlo1,dimup1
     do j=1,dimi(i,1)
        do k=1,dimi(i,2)
           write(unit,*) (item(i)%m(j,k,l), l=1,dimi(i,3))
        enddo
     enddo
     write(unit,*) ' '
  enddo
  deallocate(dimi)
  write(unit,*) '--------------------------------------'
end subroutine print_arrmat3_vec

subroutine print_arrmat2_vec(unit,item,name)
  use type_module
  use datatype
  implicit none
  integer(kind=i4_kind), intent(in)     :: unit
  character*(*),optional, intent(in)    :: name
  type(arrmat2), intent(in)             :: item(:)
  !** End of interface *****************************************
  integer(kind=i4_kind)                 :: dimlo1,dimup1,i,j,k
  integer(kind=i4_kind),allocatable     :: dimi(:,:)
  dimlo1 = lbound(item,1)
  dimup1 = ubound(item,1)
  allocate( dimi(dimup1-dimlo1+1,2))
  do i=dimlo1,dimup1
     dimi(i,1) = ubound(item(i)%m,1)
     dimi(i,2) = ubound(item(i)%m,2)
  enddo
  write(unit,*)' '
  if (present(name)) then
     write(unit,*) 'Vektor    :',name
  endif
  do i = dimlo1,dimup1
     do j=1,dimi(i,1)
        do k=1,dimi(i,2)
           write(unit,*) item(i)%m(j,k)
        enddo
     enddo
     write(unit,*) ' '
  enddo
  deallocate(dimi)
  write(unit,*) '--------------------------------------'
end subroutine print_arrmat2_vec

end module print_module

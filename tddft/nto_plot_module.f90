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
!=====================================================================
! Public interface of module
!=====================================================================
module nto_plot_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Interfere in the orbital_plot_module and convert the MO grid points
  !           into NTO grid points
  !
  !
  !  Author: Huix-Rotllant, Miquel
  !  Date: 09/2007
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
#include <def.h>
  use type_module               ! type specification parameters
  use iounitadmin_module        ! routines for I/O
  use output_module             ! contains output options
  implicit none ! by default, all the variables have to be specifically declared
  save          ! save all variables defined in this module
  private       ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------
  !------------ Declaration of constants and variables ---------------

  ! This matrix contains all the MOs that take part in the NTO
  ! expansions:
  integer(kind=i4_kind), allocatable , public :: nto_orb(:,:)
  !------------ Interface statements ---------------------------------
  !------------ public functions and subroutines ---------------------
  public nto_plot_search     ! Subroutine called by orbital_plot_module
  public nto_plot_order      ! Subroutine called by orbital_plot_module
  public nto_plot_grid       ! Subroutine called by orbital_plot_module
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !------------ Declaration of types ---------------------------------

  type, public :: coeff
     integer(i4_kind) :: index ! label of the NTO

     integer(i4_kind) :: dimocc,dimunocc ! number of occupied and
                                         ! unoccupied MO

     integer(i4_kind), allocatable :: occ_idx(:,:) ! occupied pairs
                                                   ! (index[1],symmetry[2])

     integer(i4_kind), allocatable :: unocc_idx(:,:) ! unoccupied
                                                     ! pairs
                                                     ! (index[1],symmetry[2])

     real(r8_kind), allocatable :: U(:,:) ! Expansion coefficients for
                                          ! the occupied states

     real(r8_kind), allocatable :: V(:,:) ! Expansion coefficients for
                                          ! the unoccupied states
  end type coeff

  type, public :: ntolist
     integer(i4_kind) :: irrep ! Global irrep of the NTO

     integer(i4_kind) :: spin ! Spin of the MOs in the NTO (1 for up
                              ! and 2 for down)

     integer(i4_kind) :: trans ! Spin of the transition (1 for SS and
                               ! 2 for ST)

     type(coeff), allocatable :: nto_list(:) ! NTO list of orbitals
  end type ntolist

  type(ntolist), allocatable, public :: nto_coeff(:) ! All the
                                                     ! information of
                                                     ! each NTO

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains
  !*************************************************************
  subroutine nto_plot_search(irrep,spin,list,current,num_nto,ss_or_st)
  !-------------------------------------------------------------------
    !  Purpose: Read the unformatted file "nto.dat", where the
    !           data concerning the ntos is contained. Read only the
    !           NTOs that one wants to draw
    !------------ Modules used ----------------------------------
    USE filename_module,   ONLY: outfile
    implicit none
    !------------ Declaration of formal parameters ---------------

    integer(kind=i4_kind), intent(in)  :: irrep ! Global Irrep

    integer(kind=i4_kind), intent(in)  :: spin ! Spin (not of the
                                               ! transition, but alpha
                                               ! or betha in an
                                               ! unrestricted
                                               ! calculation)

    integer(kind=i4_kind), intent(in)  :: list(:) ! List of all NTOs
                                                  ! for a the irrep

    integer(kind=i4_kind), intent(in)  :: num_nto ! index that links
                                                  ! the values for the
                                                  ! irrep in the
                                                  ! nto_coeff matrix

    integer(kind=i4_kind), intent(in)  :: current ! Total number of
                                                  ! ntos in that irrep

    integer(kind=i4_kind), intent(in)  :: ss_or_st ! Specify if one is
                                                   ! interested in
                                                   ! singlet-singlet
                                                   ! or
                                                   ! singlet-triplet

    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: io_unit
    integer(kind=i4_kind)                :: count, nto_size
    integer(kind=i4_kind)                :: dimocc,dimunocc,mindim
    integer(kind=i4_kind)                :: num_sp,num_irrep,num_ntos
    integer(kind=i4_kind)                :: a,s,idx
    integer(kind=i4_kind),allocatable    :: tmp_nto(:,:)
    integer(kind=i4_kind)                :: alloc_stat
    logical                              :: it_exists,is_opened
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    !-----------------------------!
    !Reading the data from nto.dat!
    !-----------------------------!
!
    !Determine if the file exists and it can be opened
    inquire(file=TRIM(outfile("nto.dat")), exist=it_exists)
    if(it_exists) then
        io_unit=openget_iounit(file=TRIM(outfile("nto.dat")), &
                               status='old',form='formatted')
        rewind io_unit
        inquire(io_unit,opened=is_opened)
        !Message if the file cannot be oppened
        if(.not.is_opened) call error_handler('nto_plot_module: File cannot be oppened:  nto.dat  ')
    else
        !Message if the file doesn't exist
        call error_handler('nto_plot_module: File not found:  nto.dat  ')
    end if
!
    ! Read  general  information of  the  nto.dat  file:  num_sp =  if
    ! there's SS, ST or both // num_irrep = number of irreps
    read(io_unit,*)num_sp,num_irrep
!
    !Make sure that the nto.dat file contains the necessary information
    if(num_sp.eq.2.and.ss_or_st.ne.2) then
            call error_handler('nto_plot_module: nto.dat contains no SS information')
    end if
    if(num_sp.eq.1.and.ss_or_st.ne.1) then
            call error_handler('nto_plot_module: nto.dat contains no ST information')
    end if
!
    !-------------------------------------------------------------------------------------------!
    !Jump in the nto.dat file until the information corresponding to the irrep we want is found !
    !                                                                                           !
    !  * The structure of nto.dat (if only one spin is present):                                !
    !                                                                                           !
    !                     Irrep 1                                                               !
    !                     Irrep 2                                                               !
    !                     Irrep 3                                                               !
    !                     ...                                                                   !
    !                                                                                           !
    !    - The reading is done sequencially                                                     !
    !                                                                                           !
    !  * The structure of nto.dat (if SS and ST are present):                                   !
    !                                                                                           !
    !                     Irrep 1, SS                                                           !
    !                     Irrep 1, ST                                                           !
    !                     Irrep 2, SS                                                           !
    !                     Irrep 2, ST                                                           !
    !                     Irrep 3, SS                                                           !
    !                     ...                                                                   !
    !                                                                                           !
    !    - If one is interested in SS, after reading SS, the ST must be jumped                  !
    !    - If one is interested in ST, the SS must be jumped, before ST is readed               !
    !                                                                                           !
    !-------------------------------------------------------------------------------------------!
!
    count=0
    if(num_sp.eq.1.or.num_sp.eq.2) then
       do
          count=count+1
          read(io_unit,*)num_ntos
          if(count.eq.irrep) exit
          call jump_to_next('irr',io_unit,num_ntos)
       end do
    else if(num_sp.eq.3) then
       do
          count=count+1
          read(io_unit,*)num_ntos
          if(count.eq.irrep_index(irrep,ss_or_st)) exit
          call jump_to_next('irr',io_unit,num_ntos)
       end do
    end if
!
    !Allocation of nto_coeff
    if(.not.allocated(nto_coeff)) then
          allocate(nto_coeff(num_nto),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler('nto_plot_module:  allocation failure')
    end if
    !Allocation of nto_coeff%nto_list
    if(.not.allocated(nto_coeff(current)%nto_list)) then
          allocate(nto_coeff(current)%nto_list(size(list)),stat=alloc_stat)
          if(alloc_stat/=0 ) call error_handler('nto_plot_module: allocation failure')
    end if
!
    count=0
    do idx=1,size(list)                                       !Run over all NTOs for a given irrep and spin
       nto_coeff(current)%nto_list(idx)%index=list(idx)       !Save the list of NTOs
!
       do                                                     !Search for the NTO labelled as idx in nto.dat
         count = count+1                                      !Number of NTO of the file
         if(count.gt.num_ntos) call error_handler('nto_plot_module: Some of the selected NTOs are not present in nto.dat')
         read(io_unit,*)dimocc,dimunocc                       !Number of occ and unocc MO in a NTO
     !
         if(list(idx).eq.count) then                          !If the NTO is not in the list
     !
            ! a) Contains the expansion coefficients of occupied states
            ! b) Contains the expansion coeff. for the unoccupied states
            ! c) Occ indexes and sym labels of the MOs in the expansion
            ! d) Unocc indexes and sym labels of the MOs in the expansion
            mindim=min(dimocc,dimunocc)
            allocate(nto_coeff(current)%nto_list(idx)%U(dimocc,mindim)      , & ! a)
                     nto_coeff(current)%nto_list(idx)%V(dimunocc,mindim)    , & ! b)
                     nto_coeff(current)%nto_list(idx)%occ_idx(dimocc,2)     , & ! c)
                     nto_coeff(current)%nto_list(idx)%unocc_idx(dimunocc,2) , & ! d)
                     stat=alloc_stat )
            if(alloc_stat /=0) call error_handler('nto_plot_module: allocation failure')
     !
        !------------------------!
        !For each NTO in the list!
        !------------------------!

            nto_coeff(current)%nto_list(idx)%dimocc=dimocc        !Save the number of occupied MOs
            nto_coeff(current)%nto_list(idx)%dimunocc=dimunocc    !Save the number of unoccupied MOs
            nto_coeff(current)%spin=spin                          !Save the spin of the MOs (unrestricted calculation)
            nto_coeff(current)%irrep=irrep                        !Save the irrep
            nto_coeff(current)%trans=ss_or_st                     !Save if the NTO is SS or ST
     !
            !------------------------!
            !For each occupied state !
            !------------------------!
            do a = 1,dimocc                                                           !Run over all occupied orbitals

                read(io_unit,*)nto_coeff(current)%nto_list(idx)%occ_idx(a,1)    , &   !Save the labels and expansion coefficients
                               nto_coeff(current)%nto_list(idx)%occ_idx(a,2)    , &   !for the occupied states
                              ! (nto_coeff(current)%nto_list(idx)%U(a,ap),ap=1,mindim)
                              nto_coeff(current)%nto_list(idx)%U(a,1)

         !
                if(.not.allocated(nto_orb)) then                                      !Initialize the nto_orb matrix
                    allocate(nto_orb(1,3))
                    nto_orb(1,1:2)=nto_coeff(current)%nto_list(idx)%occ_idx(a,1:2)
                    nto_orb(1,3)=spin
                    cycle
                end if
         !
                ! Save the read  MO in the nto_orb matrix  if and only
                ! if there's no other  pair (MO labe, sym label, spin)
                ! It it exists, it is not needed to put it again
                if(find_rep(nto_coeff(current)%nto_list(idx)%occ_idx(a,2),&
                            nto_coeff(current)%nto_list(idx)%occ_idx(a,1),&
                            spin)) then
                    nto_size=size(nto_orb,1)
                    allocate(tmp_nto(nto_size,3))
                    tmp_nto = nto_orb
                    deallocate(nto_orb)
                    allocate(nto_orb(nto_size+1,3))
                    nto_orb = tmp_nto
                    deallocate(tmp_nto)
                    nto_orb(nto_size+1,1:2)=nto_coeff(current)%nto_list(idx)%occ_idx(a,1:2)
                    nto_orb(nto_size+1,3)=spin
                end if
            end do
     !
            !--------------------------!
            !For each unoccupied state !
            !--------------------------!
            do s = 1,dimunocc
                read(io_unit,*)nto_coeff(current)%nto_list(idx)%unocc_idx(s,1), &      !Save the labels and expansion coefficients
                               nto_coeff(current)%nto_list(idx)%unocc_idx(s,2), &      !for the unoccupied states
                              !(nto_coeff(current)%nto_list(idx)%V(s,sp),sp=1,mindim)
                              nto_coeff(current)%nto_list(idx)%V(s,1)
         !
                ! Save the read  MO in the nto_orb matrix  if and only
                ! if there's no other  pair (MO labe, sym label, spin)
                ! It it exists, it is not needed to put it again
                if(find_rep(nto_coeff(current)%nto_list(idx)%unocc_idx(s,2),&
                            nto_coeff(current)%nto_list(idx)%unocc_idx(s,1),&
                            spin)) then
                    nto_size=size(nto_orb,1)
                    allocate(tmp_nto(nto_size,3))
                    tmp_nto = nto_orb
                    deallocate(nto_orb)
                    allocate(nto_orb(nto_size+1,3))
                    nto_orb = tmp_nto
                    deallocate(tmp_nto)
                    nto_orb(nto_size+1,1:2)=nto_coeff(current)%nto_list(idx)%unocc_idx(s,1:2)
                    nto_orb(nto_size+1,3)=spin
                end if
            end do
     !
            exit               !Once the NTO is found, exit the loop and go to the next idx index
     !
         else
            call jump_to_next('nto',io_unit,dimocc,dimunocc)        !Jump to the next NTO
         end if
       end do
    end do

    CALL returnclose_iounit(io_unit)     !Close the nto.dat file

    contains
    !*********************************************************************
      logical function find_rep(mo_irrep,mo_index,mo_spin)
        !  Purpose: Determine if (mo_irrep,mo_index,mo_spin) exist in
        !           the matrix nto_orb
        !------------ Modules used ----------------------------------
        implicit none
        !------------ Declaration of formal parameters ---------------
        integer(kind=i4_kind),intent(in)     :: mo_irrep        ! Irrep index
        integer(kind=i4_kind),intent(in)     :: mo_index        ! MO    index
        integer(kind=i4_kind),intent(in)     :: mo_spin         ! Spin  (1 for up, 2 for down)
        !** End of interface *****************************************
        !------------ Declaration of local variables -----------------
        integer(kind=i4_kind)                :: mo_idx
        logical                              :: found
        !------------ Executable code --------------------------------
        found=.true.
        do mo_idx=1,size(nto_orb,1)
            if(mo_irrep.eq.nto_orb(mo_idx,2).and.mo_spin.eq.nto_orb(mo_idx,3).and. &
               mo_index.eq.nto_orb(mo_idx,1))  found=.false.
        end do
        find_rep=found
      end function find_rep
    !*********************************************************************
      recursive subroutine jump_to_next(item,io_unit,diocc,diunocc)
        !  Purpose: Jump the io_unit file to the next irrep or nto
        !------------ Modules used ----------------------------------
        implicit none
        !------------ Declaration of formal parameters ---------------
        integer(kind=i4_kind) , intent(in)           :: io_unit        !Reading unit
        character(len=3)      , intent(in)           :: item           !What is jumped, 'irr' for irrep, 'nto' for nto
        integer(kind=i4_kind) , intent(in), optional :: diocc,diunocc  !Occupied and Unoccupied dimensions
        !** End of interface *****************************************
        !------------ Declaration of local variables -----------------
        integer(kind=i4_kind)                :: nto_j                  !loop running over ntos
        integer(kind=i4_kind)                :: dim_j                  !loop running over dimensions
        integer(kind=i4_kind)                :: dimocc,dimunocc        !Occupied and unoccupied dimensions
        !------------ Executable code --------------------------------
        select case (item)
    !
            case('irr')                                             !Jump Irrep
               do nto_j=1,diocc
                  read(io_unit,*)dimocc,dimunocc
                  call jump_to_next('nto',io_unit,dimocc,dimunocc)
               end do
    !
            case('nto')                                             !Jump nto
               do dim_j=1,diocc
                   read(io_unit,*)
               end do
               do dim_j=1,diunocc
                   read(io_unit,*)
               end do
    !
        end select
      end subroutine jump_to_next
    !*********************************************************************
      integer(kind=i4_kind) function irrep_index(irrep,ss_or_st)
        !  Purpose: Get the position of a specific irrep in the nto file
        !------------ Modules used ----------------------------------
        implicit none
        !------------ Declaration of formal parameters ---------------
        integer(kind=i4_kind) , intent(in)           :: irrep
        integer(kind=i4_kind) , intent(in)           :: ss_or_st
        !** End of interface *****************************************
        !------------ Declaration of local variables -----------------
        !------------ Executable code --------------------------------
        if(ss_or_st.eq.1) then
             irrep_index=2*irrep-1
        else
             irrep_index=2*irrep
        end if
      end function irrep_index

    end subroutine nto_plot_search
!*************************************************************
    subroutine nto_plot_grid(irrep,index,occupied,point_list,final_points)
        !  Purpose: Get the grid points for the MOs and combine them according to
        !           the NTO expansion
        !------------ Modules used ----------------------------------
        implicit none
        !------------ Declaration of formal parameters ---------------
        integer(i4_kind), intent(in)    :: irrep,index      !Irrep and Index of the NTO
        real(r8_kind)   , intent(in)    :: point_list(:,:)  !Contains the grid points of all the MOs
        logical         , intent(in)    :: occupied         !If true, calculate the grid for the NTO occupied
        real(r8_kind)   , intent(out)   :: final_points(:)  !Contains the grid points of the NTO
        !** End of interface *****************************************
        !------------ Declaration of local variables -----------------
        integer(i4_kind)                :: i_grid, i_exp
        integer(i4_kind)                :: mo_irrep,mo_index
        !------------ Executable code --------------------------------

        final_points=0
!
        if(occupied) then
           do  i_grid=1,size(point_list,1)
             do i_exp=1,nto_coeff(irrep)%nto_list(index)%dimocc
                mo_irrep = nto_coeff(irrep)%nto_list(index)%occ_idx(i_exp,2)
                mo_index = nto_coeff(irrep)%nto_list(index)%occ_idx(i_exp,1)
                final_points(i_grid)= final_points(i_grid) + nto_coeff(irrep)%nto_list(index)% &
                                      U(i_exp,1)*point_list(i_grid,get_index(mo_irrep,mo_index))
             end do
           end do
!
        else
           do  i_grid=1,size(point_list,1)
               do i_exp=1,nto_coeff(irrep)%nto_list(index)%dimunocc
                  mo_irrep = nto_coeff(irrep)%nto_list(index)%unocc_idx(i_exp,2)
                  mo_index = nto_coeff(irrep)%nto_list(index)%unocc_idx(i_exp,1)
                  final_points(i_grid)= final_points(i_grid) + nto_coeff(irrep)%nto_list(index)% &
                                        V(i_exp,1)*point_list(i_grid,get_index(mo_irrep,mo_index))
               end do
           end do
        end if
!
        contains

        integer(i4_kind) function get_index(mo_irrep,mo_label)
           !  Purpose: Find the index of the MO pair (mo_irrep, mo_label)
           !------------ Modules used ----------------------------------
           implicit none
           !------------ Declaration of formal parameters ---------------
           integer(i4_kind),intent(in) :: mo_irrep
           integer(i4_kind),intent(in) :: mo_label
           !** End of interface *****************************************
           !------------ Declaration of local variables -----------------
           integer(i4_kind)            :: count
           !------------ Executable code --------------------------------
           count = 0
           do
              if(mo_irrep.eq.nto_orb(count,2).and.mo_label.eq.nto_orb(count,1)) then
                   get_index = count
                   exit
              end if
              count = count + 1
           end do
        end function
    end subroutine nto_plot_grid
  !*************************************************************
    subroutine nto_plot_order()
        !  Purpose: Put in ascending order the MO contained in nto_orb,
        !           first by spin, second by irrep and third by label
        !------------ Modules used ----------------------------------
        implicit none
        !------------ Declaration of formal parameters ---------------
        !** End of interface *****************************************
        !------------ Declaration of local variables -----------------
        integer(i4_kind)        :: k_irrep,k_mo
        integer(i4_kind)        :: transverse(3)
!
        !Order by spin (first the ones with spin 1 and second the ones with spin 2)
!
        k_mo=1
        do k_irrep=1,size(nto_orb,1)
           if(nto_orb(k_irrep,3).eq.1) then
              transverse(:)=nto_orb(k_mo,:)
              nto_orb(k_mo,:)=nto_orb(k_irrep,:)
              nto_orb(k_irrep,:)=transverse(:)
              k_mo=k_mo+1
              cycle
           end if
        end do
!
        !Order by Symmetry label, preserving the previous order
!
        do k_irrep=1,size(nto_orb,1)
           do k_mo=k_irrep,size(nto_orb,1)
              if(nto_orb(k_mo,3).ne.nto_orb(k_irrep,3)) cycle
              if(nto_orb(k_mo,2).lt.nto_orb(k_irrep,2)) then
              transverse(:)=nto_orb(k_irrep,:)
              nto_orb(k_irrep,:)=nto_orb(k_mo,:)
              nto_orb(k_mo,:)=transverse(:)
              end if
           end do
        end do
!
        !Order by MO label, preserving the previous order
!
        do k_irrep=1,size(nto_orb,1)
             do k_mo=k_irrep,size(nto_orb,1)
                if(nto_orb(k_mo,3).ne.nto_orb(k_irrep,3)) cycle
                if(nto_orb(k_mo,2).ne.nto_orb(k_irrep,2)) cycle
                if(nto_orb(k_mo,1).lt.nto_orb(k_irrep,1)) then
                   transverse(1) = nto_orb(k_irrep,1)
                   nto_orb(k_irrep,1)=nto_orb(k_mo,1)
                   nto_orb(k_mo,1)=transverse(1)
                end if
             end do
        end do

    end subroutine nto_plot_order
  !--------------- End of module -------------------------------------
end module nto_plot_module

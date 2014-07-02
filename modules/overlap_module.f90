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
module overlap_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: 
  !        Contains the following PUBLIC variables:
  !
  !        Overlap matrix    OVERLAP
  !
  !
  !  Module called by: several
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification 
  ! Author: MM
  ! Date:   10/97
  ! Description: extension to spin orbit
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

  !------------ Modules used -----------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use datatype    ! user defined types
  use symmetry_data_module ! symmetry information

  implicit none
  private
  save
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of public constants and variables -----
  public arrmat2
  type(arrmat2), allocatable, public, protected :: overlap(:)
  type(arrmat2), allocatable, public, protected :: overlap_real(:),overlap_imag(:)

  !------------ public functions and subroutines ---------------------
  public read_overlap, dealloc_overlap


  !===================================================================
  ! End of public interface of module
  !===================================================================


  character(len=32), private :: overlap_state = 'deallocated'


!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains
  
  !*************************************************************
  subroutine read_overlap()
    !  Purpose: read in the overlap matrix
    !          Folke Noertemann   7/95
    ! subroutine called by: solve_serial,solve_para
    !------------ Modules ----------------------------------------
    use iounitadmin_module   ! provide I/O-units
    use filename_module      ! set I/O-Filenames
    use print_module         ! pretty-print interfaces
    use output_module, only:  output_overlap
    use options_module, only: options_integrals_on_file,options_spin_orbit
    use integralstore_module, only: integralstore_2cob_ol
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: i_gamma,m,n,io_u,io_u_help,i,i_meta

    ! dimensions of irreps
    ! (in order to account for SPIN ORBIT more easily)
    integer(kind=i4_kind)                :: n_irrep
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
    ! n_irrep    : number of irreps
    ! dim_irrep : number of independent functions in irrep
    !------------ Executable code ------------------------------------

    DPRINT 'ovl::read_overlap: entered, overlap_state=',overlap_state
    if( overlap_state == 'valid' )then
      WARN('no need to re-read, fix logic')
      RETURN
    endif
    ASSERT(overlap_state=='deallocated')

    ! set appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n) = ssym%dim_proj(n)
       enddo
    else
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n)  = ssym%dim(n)
       enddo
    endif

    ! sets the overlap_state to ``allocated'':
    call alloc_overlap(ssym)

    ASSERT(overlap_state=='allocated')
    ! upon exit overlap_state will be set to ``valid''

    if (options_spin_orbit) then
      ASSERT(size(overlap_real)==n_irrep)
      do i=1,n_irrep
      n = dim_irrep(i)
      ASSERT(size(overlap_real(i)%m,1)==n)
      ASSERT(size(overlap_real(i)%m,2)==n)
      ASSERT(size(overlap_imag(i)%m,1)==n)
      ASSERT(size(overlap_imag(i)%m,2)==n)
      enddo
    else
      ASSERT(size(overlap)==n_irrep)
      do i=1,n_irrep
      n = dim_irrep(i)
      ASSERT(size(overlap(i)%m,1)==n)
      ASSERT(size(overlap(i)%m,2)==n)
      enddo
    endif

    if ( output_overlap ) then
       io_u_help = openget_iounit(trim(outfile('overlapping.out')), form='formatted',status='replace')
    endif
    
    if ( options_integrals_on_file() ) then

       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          ! real part
          io_u = openget_iounit(trim(tmpfile('overlap_real.dat')), form='unformatted', status='old')
          do i_gamma = 1,n_irrep
             do m = 1,dim_irrep(i_gamma)
                do n = 1,m    
                   read(io_u) overlap_real(i_gamma)%m(m,n) 
                   overlap_real(i_gamma)%m(n,m) = overlap_real(i_gamma)%m(m,n)
                enddo
             enddo
             if ( output_overlap ) then
                call banner(io_u_help,'Overlap matrix!')
                call printout(io_u_help,overlap_real(i_gamma)%m,'Overlap_g_real')
             endif
          enddo! i_gamma loop
          call return_iounit(io_u)
          close (io_u)
          ! imaginary part
          io_u = openget_iounit(trim(tmpfile('overlap_imag.dat')), form='unformatted', status='old')
          do i_gamma = 1,n_irrep
             do m = 1,dim_irrep(i_gamma)
                do n = 1,m    
                   read(io_u) overlap_imag(i_gamma)%m(m,n) 
                   overlap_imag(i_gamma)%m(n,m) = - overlap_imag(i_gamma)%m(m,n)
                enddo
             enddo
             if ( output_overlap ) then
                call banner(io_u_help,'Overlap matrix!')
                call printout(io_u_help,overlap_imag(i_gamma)%m,'Overlap_g_imag')
             endif
          enddo! i_gamma loop
          call return_iounit(io_u)
          close (io_u)

       else ! options_spin_orbit
          io_u = openget_iounit(trim(tmpfile('overlap.dat')), form='unformatted', status='old')
          do i_gamma = 1,ssym%n_irrep
             do m = 1,ssym%dim(i_gamma)
                do n = 1,m    
                   read(io_u) overlap(i_gamma)%m(m,n) 
                   overlap(i_gamma)%m(n,m) = overlap(i_gamma)%m(m,n)
                enddo
             enddo
             if ( output_overlap ) then
                call banner(io_u_help,'Overlap matrix!')
                call printout(io_u_help,overlap(i_gamma)%m,'Overlap_g')
             endif
          enddo! i_gamma loop

          call return_iounit(io_u)
          close (io_u)
       endif ! options_spin_orbit

    else ! .not. options_integrals_on_file()

ASSERT(allocated(integralstore_2cob_ol))
       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          i_meta = 1
          do i_gamma = 1,n_irrep
             do m = 1,dim_irrep(i_gamma)
                do n = 1,m    
                   overlap_real(i_gamma)%m(m,n) = integralstore_2cob_ol(i_meta)
                   overlap_real(i_gamma)%m(n,m) = overlap_real(i_gamma)%m(m,n)
                   i_meta = i_meta + 1
                   overlap_imag(i_gamma)%m(m,n) = integralstore_2cob_ol(i_meta)
                   overlap_imag(i_gamma)%m(n,m) = - overlap_imag(i_gamma)%m(m,n)
                   i_meta = i_meta + 1
                enddo
             enddo
             if ( output_overlap ) then
                call banner(io_u_help,'Real part of Overlap matrix!')
                call printout(io_u_help,overlap_real(i_gamma)%m,'Overlap_g')
                call banner(io_u_help,'Imaginary part of Overlap matrix!')
                call printout(io_u_help,overlap_imag(i_gamma)%m,'Overlap_g')
             endif
          enddo! i_gamma loop
       else
          i_meta = 1
          do i_gamma = 1,ssym%n_irrep
             do m = 1,ssym%dim(i_gamma)
                do n = 1,m    
                   overlap(i_gamma)%m(m,n) = integralstore_2cob_ol(i_meta)
                   overlap(i_gamma)%m(n,m) = overlap(i_gamma)%m(m,n)
                   i_meta = i_meta + 1
                enddo
             enddo
             if ( output_overlap ) then
                call banner(io_u_help,'Overlap matrix!')
                call printout(io_u_help,overlap(i_gamma)%m,'Overlap_g')
             endif
          enddo! i_gamma loop
       endif ! options_spin_orbit

    endif


    if ( output_overlap ) then
       call return_iounit(io_u_help)
       close (io_u_help)
    endif

    ! deallocate appropriate dimensions of irreps
    ! (use projective irreps in case of spin orbit)
    deallocate(dim_irrep)

    ! set the global module variable:
    overlap_state ='valid'
  end subroutine read_overlap


  subroutine alloc_overlap(ssym)
    !  Purpose: allocate the appropriate space for OVERLAP
    !           FN 10/95
    use options_module, only: options_spin_orbit
    use init_module, only: init
    type(sym),intent(in)   :: ssym
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)  :: i
    !------------ Executable code ------------------------------------

    DPRINT 'ovl::alloc_overlap: entered, overlap_state=',overlap_state

    ! set the global module variable:
    overlap_state ='allocated'

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       allocate (overlap_real(ssym%n_proj_irrep),overlap_imag(ssym%n_proj_irrep))
       do i=1,ssym%n_proj_irrep
          allocate( overlap_real(i)%m(ssym%dim_proj(i),ssym%dim_proj(i)),&
               overlap_imag(i)%m(ssym%dim_proj(i),ssym%dim_proj(i)))
          call init(overlap_real(i)%m)
          call init(overlap_imag(i)%m)
       enddo
    else ! options_spin_orbit
       allocate (overlap(ssym%n_irrep))
       do i=1,ssym%n_irrep
          allocate( overlap(i)%m(ssym%dim(i),ssym%dim(i)))
       enddo
    endif
    return
  end subroutine alloc_overlap
  !*************************************************************

  !*************************************************************
  subroutine dealloc_overlap()
    implicit none
    ! *** end of interface ***

    integer(i4_kind) :: i_gamma, alloc_stat


    DPRINT 'ovl::dealloc_overlap: entered, overlap_state=',overlap_state

    ! set the global module variable:
    overlap_state ='deallocated'

    ! deallocation of overlap matrix --------------------------------
    if(allocated(overlap_real))then
       do i_gamma=1,size(overlap_real)
          deallocate(overlap_real(i_gamma)%m,overlap_imag(i_gamma)%m,STAT=alloc_stat)
          if(alloc_stat.ne.0) call error_handler &
               ("dealloc_overlap: deallocation overlap matrix failed (1)")
       enddo
       deallocate(overlap_real,overlap_imag,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("dealloc_overlap: deallocation overlap matrix failed (2)")
    endif
    if(allocated(overlap))then
       do i_gamma=1,size(overlap)
          deallocate(overlap(i_gamma)%m,STAT=alloc_stat)
          if(alloc_stat.ne.0) call error_handler &
               ("dealloc_overlap: deallocation overlap matrix failed (3)")
       enddo
       deallocate(overlap,STAT=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("dealloc_overlap: deallocation overlap matrix failed (4)")
    endif
  end subroutine dealloc_overlap
  !*************************************************************

  !--------------- End of module -------------------------------------
end module overlap_module

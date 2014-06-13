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
module int_send_aux_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
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

# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind  ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------


  !------------ Declaration of constants and variables ------------



  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: init
  public :: done
  public :: int_write_sa_int

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------


  !------------ Declaration of constants and variables ----

  logical               :: initialized = .false.
  integer(IK),parameter :: debug       = 1

  integer(IK)           :: MyIndex     = -1

  integer(IK),pointer :: offs_0(:,:),offs_L(:,:),offs_S(:,:)
  ! ...(max_bipel_index,n_irr)
  integer(IK),pointer :: dims_0(:,:),dims_L(:,:),dims_S(:,:)
  ! ...(max_bipel_index,n_irr)
  integer(IK),pointer :: irrep_size_0(:),irrep_size_L(:),irrep_size_S(:)
  ! ...(n_irr)
  integer(IK), allocatable :: qmask(:,:)
  ! ...(max_bipel,max_bipel), stores process ID which does the quadrupel


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine init()
    use error_module, only: error
    use comm_module, only: comm_myindex
    use quadrupel_fname, qfn_init=>init
    use dimensions, only:&
         & UBasDimSpor,UBasDimSporL,UBasDimSporS,&
         & n_irr=>number_of_irreps
    implicit none
    ! *** end of interface ***

    integer     :: memstat
    integer(IK) :: maxix

    call error(initialized,"int_send_aux_init: already called")

    call qfn_init()

    maxix = max_bipel_index()

    allocate(&
         & offs_0(maxix,n_irr),&
         & offs_L(maxix,n_irr),&
         & offs_S(maxix,n_irr),&
         & dims_0(maxix,n_irr),&
         & dims_L(maxix,n_irr),&
         & dims_S(maxix,n_irr),&
         & STAT=memstat)
    call error(memstat,"int_send_aux_init: alloc 1 failed")

    allocate(&
         & irrep_size_0(n_irr),&
         & irrep_size_L(n_irr),&
         & irrep_size_S(n_irr),&
         & STAT=memstat)
    call error(memstat,"int_send_aux_init: alloc 2 failed")

    call DimsAndOffs(UBasDimSpor, dims_0,offs_0,irrep_size_0)

    call DimsAndOffs(UBasDimSporL,dims_L,offs_L,irrep_size_L)

    call DimsAndOffs(UBasDimSporS,dims_S,offs_S,irrep_size_S)

    MyIndex     = comm_myindex()

    call init_qmask(maxix)

    initialized = .true.
  end subroutine init

  subroutine done()
    use comm_module
    use error_module, only: error
    use quadrupel_fname, qfn_done=>done
    implicit none
    ! *** end of interface ***

    integer :: memstat

    call error(.not.initialized,"int_send_aux_done: not initialized")

    call int_rewriteLSv2()

    !---------------------------------
    ! now actual shutdown >>>

    call qfn_done()

    deallocate(&
         & offs_0,offs_L,offs_S,&
         & dims_0,dims_L,dims_S,&
         & STAT=memstat)
    call error(memstat,"int_send_aux_done: dealloc 1 failed")

    deallocate(&
         & irrep_size_0,&
         & irrep_size_L,&
         & irrep_size_S,&
         & STAT=memstat)
    call error(memstat,"int_send_aux_done: dealloc 2 failed")

    call done_qmask()

    initialized = .false.
  end subroutine done

  subroutine init_qmask(maxix)
    use error_module
    implicit none
    integer(IK), intent(in) :: maxix
    ! *** end of interface ***
    integer(IK) :: memstat

    allocate(qmask(maxix,maxix),STAT=memstat)
    call error(memstat,"isam/init_qmask: alloc failed")
    qmask = -1
  end subroutine init_qmask

  subroutine set_qmask(ix1,ix2)
    use error_module
    implicit none
    integer(IK), intent(in) :: ix1,ix2
    ! *** end of interface ***

#ifndef FPP_NODEBUG
    if(qmask(ix1,ix2).ne.-1)then
       DPRINT MyID,'isam/set_qmask: WARNING: quadrupel(',ix1,',',ix2,') is already set'
       if(qmask(ix1,ix2).ne.MyIndex)then
          call error("isam/set_qmask: wrong process ID")
       endif
    endif
    DPRINT MyID,'isam/set_qmask: set quadrupel(',ix1,',',ix2,')=',MyIndex
#endif
    qmask(ix1,ix2) = MyIndex
  end subroutine set_qmask

  subroutine done_qmask()
    use error_module
    implicit none
    ! *** end of interface ***
    integer(IK) :: memstat

    deallocate(qmask,STAT=memstat)
    call error(memstat,"isam/done_qmask: dealloc failed")
  end subroutine done_qmask

  subroutine int_write_sa_int(sa,q,ext)
    use quadrupel_module
    use symm_adapt_int, only: sa_int_block
    use quadrupel_fname
#ifdef FPP_DEBUG
    use error_module, only: MyID
#endif
    implicit none
    type(sa_int_block), intent(in)  :: sa(:) ! (n_irr)
    type(quadrupel_type),intent(in) :: q
    character(len=3),intent(in)     :: ext
    optional :: ext
    ! *** end of interface ***

    character(len=3) :: ext_

    integer(IK) :: irr
    integer(IK) :: ix1,ix2

    DPRINT MyID,'int_write_sa_int: entered'

    ext_ = "aux"
    if(present(ext))ext_=ext

    ix1 = bipel_index(q%ua1,q%l1)
    ix2 = bipel_index(q%ua2,q%l2)

    call set_qmask(ix1,ix2)

    do irr = 1, size(sa) ! n_irr
       call int_write_sa_int_block(qfilename3("q", irr, ix1, ix2, ext_), sa(irr))
    enddo
  end subroutine int_write_sa_int

  subroutine int_write_sa_int_block(fname,b)
    ! FIXME: uses intermediate copy for reshape!
    use symm_adapt_int, only: sa_int_block
    use io, only: write_buffer
    implicit none
    character(len=*),intent(in)   :: fname
    type(sa_int_block),intent(in) :: b
    ! *** end of interface ***

    real(RK)    :: buf(size(b%re),2)

    if(size(b%re).eq.0) RETURN

    DPRINT 'int_write_sa_int_block: >',trim(fname),shape(b%re),'x2'

    buf(:,1) = reshp(b%re)
    buf(:,2) = reshp(b%im)
    call write_buffer(fname,buf)

    ! *** Internal funct ion >>> ***
  contains

    function reshp(a4d) result(a1d)
      real(RK),intent(in) :: a4d(:,:,:,:)
      real(RK)            :: a1d(size(a4d)) !<<< result
      ! *** end of interface ***

      integer(IK),parameter :: ord(4) = (/1,3,2,4/)
      integer(IK)           :: sh(4),shord(4)

      sh  = shape(a4d)

!      a1d = pack(reshape(a4d,sh(ord),ORDER=ord),.true.)
        shord = sh(ord)
      a1d = pack(reshape(a4d,shord,ORDER=ord),.true.)
    end function reshp

  end subroutine int_write_sa_int_block

  subroutine int_rewriteLSv2()
    use error_module
    use filename_module, strlen=>filename_namelengthmax
    use quadrupel_fname
    use dimensions, n_irr=>number_of_irreps
    use comm_module !, only: comm_i_am_master
    use comm, only: comm_barrier
    use msgtag_module, only: msgtag_packed_message
    use xpack
    use io, only: write_buffer, read_buffer
    implicit none
    ! *** end of interface ***

    integer          :: memstat
    integer(IK)      :: max_bipel,n_L,n_S,ix1,ix2,irr !,a1,a2,z1,z2
    real(RK),pointer :: buf(:,:,:),ints(:,:,:)
    integer(IK)      :: jx1,jx2,jrr
    logical          :: my_q
#ifndef FPP_GFORTRAN_BUGS
    real(RK) :: chksum, chksum1
#endif

    max_bipel = max_bipel_index()

    irr_: do irr=1,n_irr

       DPRINT MyID,'isam/rewriteLSv2: barrier(',irr,') ...'
       call comm_barrier()
       DPRINT MyID,'isam/rewriteLSv2:                  ... passed'

       n_L = irrep_size_L(irr)
       n_S = irrep_size_S(irr)

       if(comm_i_am_master())then
          allocate(buf(n_L,n_S,2),STAT=memstat)
          call error(memstat,"isam/rewriteLSv2: alloc buf failed")
       endif

       do ix1=1,max_bipel
          do ix2=1,max_bipel
             my_q = ( qmask(ix1,ix2).eq.MyIndex )

             if(my_q)then
                !
                ! this quadrupel is on my disk:
                call make_storage(irr,ix1,ix2)
                call read_buffer(qfilename3("q", irr, ix1, ix2, "bLS"), ints)
                !
                ! I read it in but I may be a slave, then I send it to master
                if (comm_parallel() .and. .not. comm_i_am_master()) then
                   call comm_init_send(comm_master_host,msgtag_packed_message)
                   call pck(irr)
                   call pck(ix1)
                   call pck(ix2)
                   call pck(ints)
#ifndef FPP_GFORTRAN_BUGS
                   chksum = sum(ints)
                   call pck(chksum)
#endif
                   call comm_send()
                endif
                call free_storage()
             else
                ! this is not my_quadrupel, use the iteration to get something
                ! from slaves
                ! here master receives ( something, may be not what you think it is )
                if(.not.comm_parallel())&
                     & call error("isam/rewriteLSv2: must never be here in serial run")
                if(comm_i_am_master())then
                   call comm_save_recv(comm_all_other_hosts,msgtag_packed_message)
                   call upck(jrr)
                   call upck(jx1)
                   call upck(jx2)
#ifndef FPP_NODEBUG
                   if( irr.ne.jrr .or. ix1.ne.jx1 .or. ix2.ne.jx2 )then
                      DPRINT MyID,'isam/rewriteLSv2: *****************************'
                      DPRINT MyID,'isam/rewriteLSv2: WARNING: data stream mixed up'
                      if(irr.ne.jrr)&
                           & call error(&
                           & "isam/rewriteLSv2: ERROR: irreps mixed up")
                   endif
#endif
                   call make_storage(jrr,jx1,jx2)
                   call upck(ints)

#ifndef FPP_GFORTRAN_BUGS       /* fails at -O1 with GFortran 4.6 */
                   call upck(chksum)
                   chksum1 = sum(ints) ! FIXME: force this number to 64 bits!
                   if (chksum /= chksum1) then
                      print *, 'isam/rewriteLSv2: chk1 = ', chksum, 'chk2 = ', chksum1, &
                           'diff =', chksum1 - chksum
                      ABORT("checksum failed")
                   endif
#endif
                 endif
             endif

          enddo
       enddo


       if(comm_i_am_master())then
          ! write complex matrix as type(cmatrix):
          call write_buffer(qfilename("alphap", irr, 'CMX'), buf)
          DPRINT MyID,'isam/rewriteLSv2: IRR=',irr,' SUM=',sum(buf)

          deallocate(buf,STAT=memstat)
          call error(memstat,"isam/rewriteLSv2: dealloc buf failed")
       endif
    enddo irr_

  contains

    subroutine make_storage(irr,ix1,ix2)
      implicit none
      integer(IK),intent(in) :: irr,ix1,ix2
      ! *** end of interface ***

      integer(IK) :: s1,s2,a1,a2,z1,z2

      s2 = dims_L(ix2,irr)
      a2 = offs_L(ix2,irr)
      z2 = a2 + s2  - 1

      s1 = dims_S(ix1,irr)
      a1 = offs_S(ix1,irr)
      z1 = a1 + s1 - 1

      if(comm_i_am_master())then
!!$ DPRINT MyID,'isam/rewriteLSv2: ints => buf(',a2,z2,a1,z1,'1:2)'
         ints => buf(a2:z2,a1:z1,1:2)
      else
!!$ DPRINT MyID,'isam/rewriteLSv2: allocate ints(',s2,s1,2,')'
         allocate(ints(s2,s1,2),STAT=memstat)
         call error(memstat,"isam/rewriteLSv2: alloc ints failed")
      endif
    end subroutine make_storage

    subroutine free_storage()
      implicit none
      ! *** end of interface ***
      integer(IK) :: memstat

      if (.not. comm_i_am_master()) then
         deallocate(ints,STAT=memstat)
         call error(memstat,"isam/rewriteLSv2: dealloc ints failed")
      endif
    end subroutine free_storage

  end subroutine int_rewriteLSv2

  subroutine DimsAndOffs(bd,dims,offs,irrpsize)
    use error_module
    use dimensions, only: SubSpaceDim
    use unique_atom_module, only: unique_atoms ! for checking
    use quadrupel_fname
    implicit none
    type(SubSpaceDim),intent(in) :: bd(:)
    ! ...(n_ua)
    integer(IK),intent(out)      :: dims(:,:),offs(:,:)
    ! ...(max_bipel_index,n_irr)
    integer(IK),intent(out)      :: irrpsize(:)
    ! ...(n_irr)
    ! *** end of interface ***

    integer(IK) :: n_ua,a,n_irr,irr,lmax,l,off,bipel
    integer(IK) :: max_bipel

    max_bipel = max_bipel_index()

    call error(size(offs,1)/=size(dims,1),&
         & "isam/DimsAndOffs: size(offs,2)/=size(dims,2)")
    call error(size(offs,1)/=max_bipel,&
         & "isam/DimsAndOffs: size(offs)/=max_bipel_index()")

    n_ua = size(bd)

    call error(n_ua<=0,&
         & "isam/DimsAndOffs: n_ua<=0")

    n_irr = bd(1)%n_irr

    call error(any(bd(:)%n_irr/=n_irr),&
         & "isam/DimsAndOffs: any(bd(:)%n_irr/=n_irr)")
    call error(n_irr/=size(irrpsize),&
         & "isam/DimsAndOffs: n_irr/=size(irrpsize)")
    call error(n_irr/=size(dims,2),&
         & "isam/DimsAndOffs: n_irr/=size(dims,2)")
    call error(n_irr/=size(offs,2),&
         & "isam/DimsAndOffs: n_irr/=size(offs,2)")

    do a=1,n_ua
       lmax  = bd(a)%lmax
       call error(lmax/=unique_atoms(a)%lmax_ob,&
            & "isam/DimsAndOffs: lmax/=lmax_ob")
    enddo

    do irr=1,n_irr

       off = 1
       do a=1,n_ua
          lmax  = bd(a)%lmax
          do l=0,lmax

             bipel = bipel_index(a,l)

             call error(bipel>max_bipel,&
                  & "isam/DimsAndOffs: bipel>max_bipel")

             dims(bipel,irr) = bd(a)%lm(l,irr)
             offs(bipel,irr) = off

             off = off + dims(bipel,irr)
          enddo
       enddo
       irrpsize(irr) = off - 1

       call error(irrpsize(irr)/=sum(dims(:,irr)),&
            & "isam/DimsAndOffs: irrpsize(irr)/=sum(dims(:,irr))")
    enddo

  end subroutine DimsAndOffs

  !--------------- End of module ----------------------------------

end module int_send_aux_module

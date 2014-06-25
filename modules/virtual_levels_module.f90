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
module virtual_levels_module
  !---------------------------------------------------------------
  !
  !  Purpose: Routines for sending and receiving the 
  !           eigenvectors that correspond to virtual levels
  !           This has to be done for the perturbation theory.
  !
  !  Contents:    send_eigvec_vir    sends eigenvectors of virtual
  !                                  levels to all slaves
  !
  !            receive_eigvec_vir    receives the eigenvectors of
  !                                  virtual levels from master
  !
  !
  !  Module called by: main_scf
  !
  !
  !  References: Copied from occupied_levels_module
  !
  !  Author: UB
  !  Date: 7/97
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   6/98
  ! Description: Extension to Spin Orbit
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------

#include "def.h"
  use type_module ! type specification parameters
  use datatype    ! user defined types
  use symmetry_data_module
  use eigen_data_module,only: eigvec,eigval,eigvec_real,eigvec_imag
  use occupation_module, only: n_occo
  use options_module, only: options_spin_orbit
  use comm_module
  use msgtag_module
  implicit none
  private         ! by default, all names are private
  save
  !== Interrupt end of public interface of module =================

!------------ Declaration of constants and variables ------------
  public arrmat2, arrmat3
  type(arrmat3), allocatable, public  :: eigvec_vir(:)
  type(arrmat2), allocatable, public  :: eigval_vir(:)

  ! in case of spin orbit the eigenvectors are complex
  type(arrmat2), allocatable, public  :: eigvec_vir_real(:)
  type(arrmat2), allocatable, public  :: eigvec_vir_imag(:)

  integer(kind=i4_kind) :: viralloc_stat(6)=0

!------------ public functions and subroutines ------------------
  public :: eigvec_vir_dealloc, print_viralloc
  public :: virtual_levels_bcast

!================================================================
! End of public interface of module
!================================================================

!------------ Subroutines ---------------------------------------
contains

  subroutine virtual_levels_bcast(n_vir)
    use comm_module
    implicit none
    integer(i4_kind), intent(in) :: n_vir
    ! *** end of interface ***

!!!    if(.not.comm_parallel()) RETURN
    ! FIXME: all that preparations in send/receive
    !        are they necessary in a serial run?

    ! eigvec_vir is allocated inside send_eigvec_vir only
    ! thus this call should be present if one  needs 
    ! eigvec_vir structure

    if(comm_i_am_master())then
       call send_eigvec_vir(n_vir)
    else
       call comm_save_recv(comm_master_host, msgtag_vir_levels)
       ! unpack, actualy:
       call receive_eigvec_vir()
    endif
  end subroutine virtual_levels_bcast

  subroutine send_eigvec_vir(n_vir)
    ! Purpose: send those eigenvectors which belong to the first
    !          n_vir virtual level.
    !  author: UB
    !  date  : 7/97
    !------------ Declaration of formal parameters --------------
    integer(kind=i4_kind), intent(in) :: n_vir
    !** End of interface *****************************************
    !------------ Declaration of local parameters ---------------
    integer(kind=i4_kind),allocatable :: eigvec_vir_dim(:)
    integer(kind=i4_kind)             :: i,info,is,m,mo,mv, &
                                         length,alloc_stat
    integer(kind=i4_kind)                :: n_irrep,n
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
    ! n_irrep    : number of irreps
    ! dim_irrep : number of independent functions in irrep

    external error_handler
    intrinsic min

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
  
    ! now allocate and set eigvec_vir_dim-----------------
    if(.not.allocated(eigvec_vir_dim)) &
         allocate(eigvec_vir_dim(n_irrep),stat=viralloc_stat(5))
    ASSERT(viralloc_stat(5).eq.0)
           viralloc_stat(5)=1
    ! find out how many virtual levels per IRREP
    ! we have (spin avaraged) and pack them into the buffer
    if (ssym%n_spin>1) then
       eigvec_vir_dim = dim_irrep-min(n_occo(1,:),n_occo(2,:))
    else
       eigvec_vir_dim = dim_irrep-n_occo(1,:)
    endif

    eigvec_vir_dim = min(n_vir, eigvec_vir_dim)
    !-----------------------------------------------------

    ! re-allocate eigvec_vir and eigval_vir -----------------
    ! because the number of virtual orbital may have changed!
    ! first make sure that they have been deallocated before
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       if(allocated(eigvec_vir_real)) then 
          do i=1,n_irrep
             deallocate(eigvec_vir_real(i)%m,eigvec_vir_imag(i)%m)
          end do
          deallocate(eigvec_vir_real,eigvec_vir_imag,stat=alloc_stat)
          if(alloc_stat.ne.0) call error_handler &
               ("SEND_EIGVEC_vir: deallocation 1 failed")
       end if
    else
       if(allocated(eigvec_vir)) then 
          do i=1,n_irrep
             deallocate(eigvec_vir(i)%m,stat=viralloc_stat(2))
             ASSERT(viralloc_stat(2).eq.0)
          end do
          deallocate(eigvec_vir,stat=viralloc_stat(1)) !!in send_eigvec_vir
          ASSERT(viralloc_stat(1).eq.0)
       end if
    endif
    if(allocated(eigval_vir)) then 
       do i=1,size(eigval_vir)
          deallocate(eigval_vir(i)%m,stat=viralloc_stat(4)) !!in send_eigvec_vir
          ASSERT(viralloc_stat(4).eq.0)
       end do
       deallocate(eigval_vir,stat=viralloc_stat(3))
          ASSERT(viralloc_stat(3).eq.0)
    end if
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       if(.not.allocated(eigvec_vir_real)) then
          allocate(eigvec_vir_real(n_irrep),eigvec_vir_imag(n_irrep),stat=alloc_stat)
          if(alloc_stat.ne.0) call error_handler &
               ("SEND_EIGVEC_vir: allocation 4 failed")         
          do i=1,n_irrep
             allocate(eigvec_vir_real(i)%m(dim_irrep(i) &
                  ,eigvec_vir_dim(i))  &
                  ,eigvec_vir_imag(i)%m(dim_irrep(i) &
                  ,eigvec_vir_dim(i))  &
                  ,stat=alloc_stat)
             if(alloc_stat.ne.0) call error_handler &
                  ("SEND_EIGVEC_vir: allocation 5 failed")            
             eigvec_vir_real(i)%m = 0.0_r8_kind
             eigvec_vir_imag(i)%m = 0.0_r8_kind
          enddo
       endif
    else ! options_spin_orbit
       if(.not.allocated(eigvec_vir)) then
          allocate(eigvec_vir(n_irrep),stat=viralloc_stat(1))
          ASSERT(viralloc_stat(1).eq.0)
                 viralloc_stat(1)=1
          do i=1,n_irrep
             allocate(eigvec_vir(i)%m(dim_irrep(i),eigvec_vir_dim(i),ssym%n_spin)  &
                  ,stat=viralloc_stat(2))
             ASSERT(viralloc_stat(2).eq.0)
                    viralloc_stat(2)=1
             eigvec_vir(i)%m = 0.0_r8_kind
          enddo
       endif
    endif! options_spin_orbit

    if(.not.allocated(eigval_vir)) then
       allocate(eigval_vir(n_irrep),stat=viralloc_stat(3))
       ASSERT(viralloc_stat(3).eq.0)
              viralloc_stat(3)=1
       do i=1,n_irrep
          allocate(eigval_vir(i)%m(eigvec_vir_dim(i) &
               ,ssym%n_spin),stat=viralloc_stat(4))
          ASSERT(viralloc_stat(4).eq.0)
                 viralloc_stat(4)=1
          eigval_vir(i)%m = 0.0_r8_kind
       enddo
    endif
    !-------------------------------------------------------
    
    ! now copy those eigenvectors into eigvec_vir 
    ! which belong to an virtual level
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do i=1,n_irrep
          m = eigvec_vir_dim(i)
          mo = n_occo(1,i)
          mv = MIN(m,dim_irrep(i)-mo)
          eigvec_vir_real(i)%m(:,1:mv) = eigvec_real(i)%m(:,mo+1:mo+mv)
          eigvec_vir_imag(i)%m(:,1:mv) = eigvec_imag(i)%m(:,mo+1:mo+mv)
          eigval_vir(i)%m(1:mv,1) = eigval(i)%m(mo+1:mo+mv,1)
       enddo
    else
       do i=1,n_irrep
          m = eigvec_vir_dim(i)
          do is=1,ssym%n_spin
             mo = n_occo(is,i)
             mv = MIN(m,dim_irrep(i)-mo)
             eigvec_vir(i)%m(:,1:mv,is) = eigvec(i)%m(:,mo+1:mo+mv,is)
             eigval_vir(i)%m(1:mv,is) = eigval(i)%m(mo+1:mo+mv,is)
          enddo
       enddo
    endif
    
    if (comm_parallel()) then ! parallel-loop --------------
       ! start packing
       call comm_init_send(comm_all_other_hosts, msgtag_vir_levels)
       call commpack(eigvec_vir_dim(1),n_irrep,1,info)
       if(info.ne.0) call error_handler &
            ("SEND_EIGVEC_vir: packing failed")

       if (options_spin_orbit) then
          !
          ! SPIN ORBIT
          !
          ! now pack the variable eigvec_vir into the comm buffer
          do i=1,n_irrep
             length=eigvec_vir_dim(i)*dim_irrep(i)*ssym%n_spin
             if (length.gt.0) then
                call commpack(eigvec_vir_real(i)%m(1,1),length,1,info)
                if(info.ne.0) call error_handler &
                     ("SEND_EIGVEC_vir: packing failed")
                call commpack(eigvec_vir_imag(i)%m(1,1),length,1,info)
                if(info.ne.0) call error_handler &
                     ("SEND_EIGVEC_vir: packing failed")
             endif
          enddo
       else
          ! now pack the variable eigvec_vir into the comm buffer
          do i=1,n_irrep
             length=eigvec_vir_dim(i)*dim_irrep(i)*ssym%n_spin
             if (length > 0) then
                call commpack(eigvec_vir(i)%m(1,1,1),length,1,info)
                if(info.ne.0) call error_handler &
                     ("SEND_EIGVEC_vir: packing failed")
             endif
          enddo
       endif

       ! now pack the variable eigval_vir into the comm buffer
       do i=1,n_irrep
          length=eigvec_vir_dim(i)*ssym%n_spin
          if (length > 0) then
             call commpack(eigval_vir(i)%m(1,1),length,1,info)
          ASSERT(info.eq.0)
          endif
       enddo

       call comm_send()
       
    endif! End parallel loop ---------------------------------

    deallocate(eigvec_vir_dim,stat=viralloc_stat(5))
    ASSERT(viralloc_stat(5).eq.0)

    deallocate(dim_irrep,stat=viralloc_stat(6))
    ASSERT(viralloc_stat(6).eq.0)
  end subroutine send_eigvec_vir

  
  subroutine receive_eigvec_vir
    !** End of interface *****************************************
    integer(kind=i4_kind) :: i,n,n_irrep,ispin,alloc_stat,info,length
    integer(kind=i4_kind),allocatable :: eigvec_vir_dim(:)

    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
    ! n_irrep    : number of irreps
    ! dim_irrep : number of independent functions in irrep

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep=ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n) = ssym%dim_proj(n)
       enddo
    else
       n_irrep=ssym%n_irrep
       allocate(dim_irrep(n_irrep))
       do n=1,n_irrep
          dim_irrep(n)  = ssym%dim(n)
       enddo
    endif
    ispin=ssym%n_spin

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       if(.not.allocated(eigvec_vir_real)) then
          allocate(eigvec_vir_real(n_irrep),eigvec_vir_imag(n_irrep),eigval_vir(n_irrep), &
               stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_vir: allocation (1) failed')
       else
          do i=1,n_irrep
             deallocate(eigvec_vir_real(i)%m,eigvec_vir_imag(i)%m,eigval_vir(i)%m)
          end do
       end if
    else ! options_spin_orbit 
       if(.not.allocated(eigvec_vir)) then
          allocate(eigvec_vir(n_irrep),eigval_vir(n_irrep), &
               stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_vir: allocation (1) failed')
       else
          do i=1,n_irrep
             deallocate(eigvec_vir(i)%m,eigval_vir(i)%m)
          end do
       end if
    endif ! options_spin_orbit 

    allocate(eigvec_vir_dim(n_irrep),stat=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ('receive_eigvec_vir: allocation (2'') failed')
    call communpack(eigvec_vir_dim(1),n_irrep,1,info)
    if(info/=0) then
       call error_handler('receive_eigvec_vir: unpacking failed')
    endif

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do i=1,n_irrep
          allocate(eigvec_vir_real(i)%m(dim_irrep(i),eigvec_vir_dim(i)),&
               eigvec_vir_imag(i)%m(dim_irrep(i),eigvec_vir_dim(i)) &
               ,stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_vir: allocation 3 failed')
          allocate(eigval_vir(i)%m(eigvec_vir_dim(i),ispin)&
               ,stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_vir: allocation 4 failed')
       enddo

       do i=1,n_irrep
          length=eigvec_vir_dim(i)*dim_irrep(i)*ispin
          if(length/=0) then
             call communpack(eigvec_vir_real(i)%m(1,1),length,1,info)
             call communpack(eigvec_vir_imag(i)%m(1,1),length,1,info)
          end if
       enddo
    else ! options_spin_orbit
       do i=1,n_irrep
          allocate(eigvec_vir(i)%m(dim_irrep(i),eigvec_vir_dim(i),ispin)&
               ,stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_vir: allocation 3 failed')
          allocate(eigval_vir(i)%m(eigvec_vir_dim(i),ispin)&
               ,stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_vir: allocation 4 failed')
       enddo

       do i=1,n_irrep
          length=eigvec_vir_dim(i)*dim_irrep(i)*ispin
          if (length > 0) then
             call communpack(eigvec_vir(i)%m(1,1,1),length,1,info)
          end if
       enddo
    endif! options_spin_orbit

    do i=1,n_irrep
       length=eigvec_vir_dim(i)*ispin
       if (length > 0) then
          call communpack(eigval_vir(i)%m(1,1),length,1,info)
       end if
    enddo

    deallocate( eigvec_vir_dim )

    ! deallocate dimensions of irreps
    deallocate(dim_irrep)
  end subroutine receive_eigvec_vir

  !***************************************************************
    
  subroutine eigvec_vir_dealloc(context)
    ! Purpose: deallocates the following variables:
    !          -eigvec_vir
    !          -eigval_vir
    ! Subroutine called by:  xc_hamiltonian,main_slave,main_gradient
    use interfaces, only: IPARA, IMAST
    implicit none
    integer(i4_kind), intent(in) :: context
    !** End of interface *****************************************
    ! ---------- declaration of local variables --------------
    integer(kind=i4_kind)    :: alloc_stat,i
    ! --------- executable code ------------------------------

    if(IAND(context,IPARA)==IMAST)then
       ! tell slaves to call me
       DPRINT ' eigvec_vir_dealloc comm_init_send + comm_send'
       if (comm_parallel() .and. comm_i_am_master()) then
          call comm_init_send(comm_all_other_hosts,msgtag_eigvec_vir_dealloc)
          call comm_send()
       endif
       DPRINT 'done'
    endif

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       do i=1,size(eigvec_vir_real) ! n_irrep
          deallocate(eigvec_vir_real(i)%m,eigvec_vir_imag(i)%m,STAT=alloc_stat)
          if (alloc_stat /= 0 ) call error_handler &
               ("eigvec_vir_dealloc : deallocation (1a) failed")
       enddo
       deallocate(eigvec_vir_real,eigvec_vir_imag,STAT=alloc_stat)
       if (alloc_stat /= 0 ) call error_handler &
            ("eigvec_vir_dealloc : deallocation (1) failed")
    else
     if(allocated(eigvec_vir)) then
       do i=1,size(eigvec_vir) ! n_irrep
          deallocate(eigvec_vir(i)%m,STAT=viralloc_stat(2))
          ASSERT(viralloc_stat(2).eq.0)
       enddo
       deallocate(eigvec_vir,STAT=viralloc_stat(1))  !! in eigvec_vir_dealloc
       ASSERT(viralloc_stat(1).eq.0)
     endif
    endif

    ! SO and NO-SO:
    if( allocated(eigval_vir) )then
      do i=1,size(eigval_vir)
        deallocate(eigval_vir(i)%m,STAT=viralloc_stat(4))
        ASSERT(viralloc_stat(4).eq.0)
      enddo
      deallocate(eigval_vir,STAT=viralloc_stat(3))
      ASSERT(viralloc_stat(3).eq.0)
    endif
  end subroutine eigvec_vir_dealloc

  subroutine print_viralloc()
    !
    ! Debug prints
    !
    implicit none
    ! *** end of interface ***

    integer(kind=i4_kind) :: i

    do i = 1, size(viralloc_stat)
        if(viralloc_stat(i).ne.0) print*,'print_viralloc: viralloc_stat ne 0', i
    enddo
  end subroutine print_viralloc
    
!--------------- End of module ----------------------------------
 end module virtual_levels_module

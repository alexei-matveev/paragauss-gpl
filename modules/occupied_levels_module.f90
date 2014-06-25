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
module occupied_levels_module
  !---------------------------------------------------------------
  !
  !  Purpose: Routines for sending and receiving the
  !           eigenvectors that correspond to occupied levels
  !           This has to be done to compute the density on
  !           the grid and the density matrix on each slave.
  !           Currently the first SCF-Cycle skips this step
  !           leaving out the XC-hamiltonian completely.
  !           Lateron it is planned to either start with
  !           a set of old eigenvectors (geometry optimization)
  !           or an empty set.
  !
  !  Contents:    send_eigvec_occ    sends eigenvectors of occupied
  !                                  levels to all slaves
  !
  !            receive_eigvec_occ    receives the eigenvectors of
  !                                  occupied levels from master
  !
  !            update_eigvec_occ     updates occupied eigenvectors
  !                                  changed by perturbation theory
  !
  !  Module called by: main_scf
  !
  !
  !  References: ...
  !
  !  Author: MS,FN
  !  Date: 1/96
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: MS
  ! Date:   3/97
  ! Description: If subroutines are called with optional argument eigval_also,
  !              not only eigenvectors, but also eigenvalues will be copied
  !              to slaves (the same holds for the perturbation theory, UB)
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   7/97
  ! Description: Pertrubation theory on the slaves is based on the data
  !              in the occupied_levels_module and the virtual_levels_module
  !              with the latter only holding those virtual states necessary
  !              for the perturbation theory. For that purpose it is important
  !              that the index of the occupied states in not altered by
  !              empty states below the highest (partially) occupied orbital
  !              of each IRREP and spin.
  !              New routines for updating the occupied eignefunctions which
  !              are rotated by the perturbation theory.
  !
  ! Modification (Please copy before editing)
  ! Author: MM
  ! Date:   6/98
  ! Description: Extension to Spin Orbit
  !
  ! Modification send and receive routines -> sndrcv
  ! Author: AN
  ! Date:   8/10
  ! Description: summarizing the seperated send and receive routines in one
  !              routine for parallel context, add also a wrapper to enter
  !              this routine from a master only context
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------

# include "def.h"
  use type_module ! type specification parameters
  use datatype    ! user defined types
  use symmetry_data_module
  implicit none
  private         ! by default, all names are private
  save
  !== Interrupt end of public interface of module =================

!------------ Declaration of constants and variables ------------
  public arrmat2,arrmat3
  type(arrmat3), allocatable, public, target, protected :: eigvec_occ(:)
  ! eigvec_occ(i_irrep)%m(i_basis,i_occ,i_spin)
  type(arrmat2), allocatable, target, public, protected :: eigval_occ(:)
  ! eigval_occ(i_irrep)%m(i_occ,i_spin)
  type(arrmat2), allocatable, target, public, protected :: occ_num_occ(:)
  ! occ_num_occ(i_irrep)%m(i_occ,i_spin)
  type(arrmat2), allocatable, public, protected :: occ_num_core_occ(:)
  ! occ_num_core_occ(i_irrep)%m(i_occ,i_spin)

  ! in the case of spin orbit eigenvectors are complex
  type(arrmat2), allocatable, target, public, protected :: eigvec_occ_real(:)
  type(arrmat2), allocatable, target, public, protected :: eigvec_occ_imag(:)

  !
  ! FIXME: target is here only to make things like
  !
  !     v => eigvec_occ_real(irrep)%m
  !     n => occ_num_occ(irrep)%m
  !
  ! compile also when arrmat2 is declared with
  ! allocatable component instead of pointer.
  !

!------------ public functions and subroutines ------------------
  public :: eigvec_occ_dealloc!()
  public sndrcv_eigvec_occ, sndrcv_eigvec_occ1

  public :: update_eigvec_occ!(), to be executed in parallel context

!================================================================
! End of public interface of module
!================================================================

!------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine sndrcv_eigvec_occ1()
    !  Purpose: sometimes eigvec are sended from a complete master only
    !           context, this functions enables the usage of sndrcv in
    !             it
    !             master sends command to join in his subroutine to slaves
    !             slaves will also enter this routine
    !------------ Modules used ------------------- ---------------
    use comm_module, only: comm_init_send, comm_send, comm_i_am_master, &
        comm_all_other_hosts
    use msgtag_module, only: msgtag_occ_levels
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code --------------------------------

    if (comm_i_am_master()) then
      call comm_init_send(comm_all_other_hosts, msgtag_occ_levels)
      call comm_send()
    end if
    call sndrcv_eigvec_occ()
  end subroutine sndrcv_eigvec_occ1

  !*************************************************************
  subroutine sndrcv_eigvec_occ()
    !  Purpose: summarizes sending and receiving of the eigenvectors
    !           of the occupied levels
    !------------ Modules used ------------------- ---------------
    use comm
    use operations_module, only: operations_core_density
    use options_module, only: options_spin_orbit
    use eigen_data_module, only: eigvec, eigval, eigvec_real, eigvec_imag
    use occupation_module, only: occ_num, n_occo, alloc_n_occo, &
        occ_num_core, n_occo_core
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind) :: i, m
    integer(i4_kind) :: rank
    !------------ Executable code --------------------------------

    rank = comm_rank()

    ! only for slaves, master is assumed to have that:
    if(.not. allocated(n_occo))then
      ASSERT(rank > 0)
      call alloc_n_occo(ssym)
    endif

    if( comm_parallel() )then ! parallel loop for giving n_occo to slaves
      call comm_bcast(n_occo)
      if (operations_core_density) then
        call comm_bcast(n_occo_core)
      endif
    endif

    ! Now slaves and master have everything for allocation
    call eigvec_occ_alloc(n_occo, n_occo_core, ssym)

    if( rank == 0 )then ! master build up eigvec occ
      ! now copy those eigenvectors into eigvec_occ
      ! which belong to an occupied level
      ! do the same with occ_num.

      ! copy eigenvalues and occupation numbers in any case:
      do i=1, size(eigval_occ)
         ! number of occupied orbitals:
         m = size(eigval_occ(i)%m, 1)

         ! second dimension is spin, either 1:1 or 1:2:
         eigval_occ(i)%m(:, :) = eigval(i)%m(1:m, :)
         occ_num_occ(i)%m(:, :) = occ_num(i)%m(1:m, :)
      enddo

      ! copy eigenvectors:
      if( options_spin_orbit )then
         !
         ! SPIN ORBIT
         !
         do i=1, size(eigvec_occ_real)
            ! number of occupied orbitals:
            m = size(eigvec_occ_real(i)%m, 2)

            eigvec_occ_real(i)%m(:, :) = eigvec_real(i)%m(:, 1:m)
            eigvec_occ_imag(i)%m(:, :) = eigvec_imag(i)%m(:, 1:m)
         enddo
      else
         do i=1, size(eigvec_occ)
            ! number of occupied orbitals:
            m = size(eigvec_occ(i)%m, 2)

            ! third dimension is spin, either 1:1 or 1:2:
            eigvec_occ(i)%m(:, :, :) = eigvec(i)%m(:, 1:m, :)
         enddo
      endif

      if( operations_core_density )then
        do i = 1, size(occ_num_core_occ)
          ! number of occupied orbitals:
          m = size(occ_num_core_occ(i)%m, 1)

          ! second dimension is spin, either 1:1 or 1:2:
          occ_num_core_occ(i)%m(:, :) = occ_num_core(i)%m(1:m, :)
        enddo
      endif
    endif ! master build eigvec occ

    ! if not parallel, then there is no need for the send of the slaves
    ! thus return point serial run
    if( .not. comm_parallel() ) return

    ! first the eigvectors and eigvalues (complex ones if spin orbits)
    do i=1, size(eigval_occ)
      call comm_bcast(eigval_occ(i)%m)
      if( options_spin_orbit )then
        call comm_bcast(eigvec_occ_real(i)%m)
        call comm_bcast(eigvec_occ_imag(i)%m)
      else
        call comm_bcast(eigvec_occ(i)%m)
      endif
    enddo

    ! next the occupation numbers
    do i=1, size(occ_num_occ)
      call comm_bcast(occ_num_occ(i)%m)
    enddo

    ! is this still used somewhere?
    if( operations_core_density )then
      do i=1, size(occ_num_core_occ)
        call comm_bcast(occ_num_core_occ(i)%m)
      enddo
    endif
  end subroutine sndrcv_eigvec_occ

  ! ***************************************************************

  subroutine eigvec_occ_alloc(n_occo, n_occo_core, ssym)
    !  Purpose: Does all the (re-)allocation needed for the occupied
    !           orbitals, this should be done on master and slaves
    !------------ Modules used ------------------- ---------------
    use operations_module, only: operations_core_density
    use options_module, only: options_spin_orbit
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent( in  ) :: n_occo(:,:), n_occo_core(:,:)
    type(sym),             intent( in  ) :: ssym
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)                :: alloc_stat
    integer(kind=i4_kind)                :: n_irrep
    integer(kind=i4_kind), allocatable   :: dim_irrep(:)
    integer(kind=i4_kind)                :: i, n
    integer(kind=i4_kind), allocatable   :: eigvec_occ_dim(:)
    !------------ Executable code --------------------------------
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       n_irrep = ssym%n_proj_irrep
       allocate(dim_irrep(n_irrep), stat=alloc_stat)
           ASSERT(alloc_stat==0)
       do n=1,n_irrep
          dim_irrep(n) = ssym%dim_proj(n)
       enddo
    else
       n_irrep = ssym%n_irrep
       allocate(dim_irrep(n_irrep), stat=alloc_stat)
           ASSERT(alloc_stat==0)
       do n=1,n_irrep
          dim_irrep(n)  = ssym%dim(n)
       enddo
    endif

    ! now allocate and set eigvec_occ_dim-----------------
    allocate(eigvec_occ_dim(n_irrep),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    ! find out how many occupied levels per IRREP
    ! we have (spin avaraged) and pack them into the buffer
    if (ssym%n_spin>1) then
       if (operations_core_density) then
          eigvec_occ_dim = max(n_occo(1,:),n_occo(2,:), &
                               n_occo_core(1,:),n_occo_core(2,:))
       else
          eigvec_occ_dim = max(n_occo(1,:),n_occo(2,:))
       endif
    else
       if (operations_core_density) then
          eigvec_occ_dim = max(n_occo(1,:),n_occo_core(1,:))
       else
          eigvec_occ_dim = n_occo(1,:)
       endif
    endif

    ! re-allocate eigvec_occ, eigval_occ and occ_num_occ -----------------
    ! because the number of occupied orbital may have changed!
    ! first make sure that they have been deallocated before
    call eigvec_occ_dealloc()

    ! now allocate
    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       allocate(eigvec_occ_real(n_irrep),&
            eigvec_occ_imag(n_irrep), stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_irrep
          allocate(eigvec_occ_real(i)%m(dim_irrep(i) &
               ,eigvec_occ_dim(i))  &
               ,stat=alloc_stat)
          ASSERT(alloc_stat==0)
          allocate(eigvec_occ_imag(i)%m(dim_irrep(i) &
               ,eigvec_occ_dim(i))  &
               ,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
    else ! options_spin_orbit
       allocate(eigvec_occ(n_irrep), stat=alloc_stat)
       ASSERT(alloc_stat==0)
       do i=1,n_irrep
          allocate(eigvec_occ(i)%m(dim_irrep(i) &
               ,eigvec_occ_dim(i),ssym%n_spin)  &
               ,stat=alloc_stat)
          ASSERT(alloc_stat==0)
       enddo
    endif ! options_spin_orbit

    allocate(occ_num_occ(n_irrep), stat=alloc_stat)
    ASSERT(alloc_stat==0)
    do i=1,n_irrep
       allocate(occ_num_occ(i)%m(eigvec_occ_dim(i) &
            ,ssym%n_spin), stat=alloc_stat)
       ASSERT(alloc_stat==0)
    enddo

    if (operations_core_density) then
      allocate(occ_num_core_occ(n_irrep))
      ASSERT(alloc_stat==0)
      do i=1,n_irrep
         allocate(occ_num_core_occ(i)%m(eigvec_occ_dim(i) &
              ,ssym%n_spin), stat=alloc_stat)
         ASSERT(alloc_stat==0)
      enddo
    endif

    allocate(eigval_occ(n_irrep), stat=alloc_stat)
    ASSERT(alloc_stat==0)
    do i=1,n_irrep
       allocate(eigval_occ(i)%m(eigvec_occ_dim(i) &
            ,ssym%n_spin), stat=alloc_stat)
       ASSERT(alloc_stat==0)
    enddo

    deallocate(eigvec_occ_dim, stat=alloc_stat)
    ASSERT(alloc_stat==0)
  end subroutine eigvec_occ_alloc

    subroutine eigvec_occ_dealloc()
      ! Purpose: deallocates the following variables:
      !          -eigvec_occ
      !          -occ_num_occ
      !          -eigval_occ
      !          -occ_num_core_occ
      !      if they are allocated at the moment
      ! Subroutine called by:  xc_hamiltonian,main_slave,main_gradient
      implicit none
      ! ---------- declaration of formal parameters ------------
      !** End of interface *****************************************
      ! ---------- declaration of local variables --------------
      integer(kind=i4_kind)    :: alloc_stat,i
      ! --------- executable code ------------------------------

      if(allocated(eigvec_occ)) then
          do i=1, size(eigvec_occ)
             deallocate(eigvec_occ(i)%m, STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          enddo
          deallocate(eigvec_occ, STAT=alloc_stat)
          ASSERT(alloc_stat==0)
      endif

      if(allocated(eigvec_occ_real)) then
          do i=1, size(eigvec_occ_real)
             deallocate(eigvec_occ_real(i)%m, eigvec_occ_imag(i)%m, STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          enddo
          deallocate(eigvec_occ_real, eigvec_occ_imag, STAT=alloc_stat)
          ASSERT(alloc_stat==0)
      endif

      if(allocated(eigval_occ)) then
          do i=1, size(eigval_occ)
             deallocate(eigval_occ(i)%m, STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          enddo
          deallocate(eigval_occ, STAT=alloc_stat)
          ASSERT(alloc_stat==0)
      endif

      if(allocated(occ_num_occ)) then
          do i=1, size(occ_num_occ)
             deallocate(occ_num_occ(i)%m, STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          enddo
          deallocate(occ_num_occ, STAT=alloc_stat)
          ASSERT(alloc_stat==0)
      endif

      if(allocated(occ_num_core_occ)) then
          do i=1, size(occ_num_core_occ)
             deallocate(occ_num_core_occ(i)%m, STAT=alloc_stat)
             ASSERT(alloc_stat==0)
          enddo
          deallocate(occ_num_core_occ, STAT=alloc_stat)
          ASSERT(alloc_stat==0)
      endif
    end subroutine eigvec_occ_dealloc

  !*************************************************************

  subroutine update_eigvec_occ()
    !
    ! Modifies eigenvectors in "eigvec_occ", I assume.
    !
    ! It was a clever trick to make it look like a plain send/recv
    ! combination.
    !
    ! To be called from a parallel context.
    !
    use comm, only: comm_rank
    use comm_module, only: comm_save_recv, comm_master_host
    use msgtag_module, only: msgtag_rot_levels
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        !
        ! FIXME: its effect in serial runs is non-trivial:
        !
        call send_eigvec_rot()
    else
        !
        ! FIXME: slaves wait unconditionally, so send above should
        !        be unconditional too:
        !
        call comm_save_recv(comm_master_host, msgtag_rot_levels)
        call receive_eigvec_rot()
    endif
  end subroutine update_eigvec_occ

  subroutine send_eigvec_rot
    ! Purpose: send those rotated eigenvectors which belong to an
    !          occupied level.
    !
    ! FIXME: This may do more than just a send. I am afraid its
    !        effect is non-trivial in serial run too.
    !
    use comm_module
    use msgtag_module, only: msgtag_rot_levels
    use options_module, only: options_spin_orbit
    use eigen_data_module, only: eigvec, eigvec_real, eigvec_imag
    !** End of interface *****************************************
    !
    !  author: UB (based on send_eigvec_occ of FN)
    !  date  : 7/97
    ! --- Modules ---
    use options_module, only: options_xcmode, &
                              xcmode_exchange_fit, xcmode_numeric_exch
    use occupation_module, only: n_occo, n_rot
    !------------ Declaration of local parameters ---------------
    logical                           :: send
    integer(kind=i4_kind)             :: i,info,is,n,length,alloc_stat
    real(kind=r8_kind), allocatable   :: eigvec_rot(:,:,:)
    ! in case of Spin Orbit eigvec_rot is complex
    type(arrmat2),allocatable         :: eigvec_rot_real(:),eigvec_rot_imag(:)
    external error_handler

    !
    ! Historically send was conditional in regular branch
    ! like this (that is almost always performed):
    !
!   send = comm_parallel() .and. &
!        ( options_xcmode() == xcmode_exchange_fit .or. &
!          options_xcmode() == xcmode_numeric_exch )
    send = comm_parallel()

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       ! allocate eigvec_rot ---------------------------------
       if (comm_parallel()) then
          allocate(eigvec_rot_real(ssym%n_proj_irrep),&
               eigvec_rot_imag(ssym%n_proj_irrep),stat=alloc_stat)
          if(alloc_stat.ne.0) call error_handler &
               ("SEND_EIGVEC_ROT: allocation 1 failed")
          do i=1,ssym%n_proj_irrep
             allocate(eigvec_rot_real(i)%m(ssym%dim_proj(i),2),eigvec_rot_imag(i)%m(ssym%dim_proj(i),2),stat=alloc_stat)
             if(alloc_stat.ne.0) call error_handler &
                  ("SEND_EIGVEC_ROT: allocation 2 failed")
             eigvec_rot_real(i)%m = 0.0_r8_kind
             eigvec_rot_imag(i)%m = 0.0_r8_kind
          enddo
       endif

       ! now copy those eigenvectors into eigvec_occ and eigvec_rot
       ! which have been changed by the perturbation theory
       do i=1,ssym%n_proj_irrep
          n=n_rot(i,1)
          if (n>0.and.n<=n_occo(1,i)) then
             eigvec_occ_real(i)%m(:,n) = eigvec_real(i)%m(:,n)
             eigvec_occ_imag(i)%m(:,n) = eigvec_imag(i)%m(:,n)
             if (n<n_occo(1,i)) then
                eigvec_occ_real(i)%m(:,n+1) = eigvec_real(i)%m(:,n+1)
                eigvec_occ_imag(i)%m(:,n+1) = eigvec_imag(i)%m(:,n+1)
             endif
             if (comm_parallel()) then
                eigvec_rot_real(i)%m(:,1) = eigvec_real(i)%m(:,n)
                eigvec_rot_imag(i)%m(:,1) = eigvec_imag(i)%m(:,n)
                if (n<n_occo(1,i)) then
                   eigvec_rot_real(i)%m(:,2) = eigvec_real(i)%m(:,n+1)
                   eigvec_rot_imag(i)%m(:,2) = eigvec_imag(i)%m(:,n+1)
                endif
             endif
          endif
       enddo

       if (comm_parallel()) then ! parallel-loop --------------
          ! start packing
          call comm_init_send(comm_all_other_hosts,msgtag_rot_levels)
          call commpack(n_rot(1,1),ssym%n_proj_irrep*ssym%n_spin,1,info)
          if(info.ne.0) call error_handler &
               ("SEND_EIGVEC_ROT: packing failed")
          do i=1,ssym%n_proj_irrep
             length=ssym%dim_proj(i)*2*ssym%n_spin
             if (length > 0) then
                call commpack(eigvec_rot_real(i)%m(1,1),length,1,info)
                if(info.ne.0) call error_handler &
                     ("SEND_EIGVEC_ROT: packing failed")
                call commpack(eigvec_rot_imag(i)%m(1,1),length,1,info)
                if(info.ne.0) call error_handler &
                     ("SEND_EIGVEC_ROT: packing failed")
             endif
          enddo

          call comm_send()

       endif! End parallel loop ---------------------------------

       ! Now deallocate eigvec_rot again
       if (comm_parallel()) then
          do i=1,ssym%n_proj_irrep
             deallocate(eigvec_rot_real(i)%m,eigvec_rot_imag(i)%m,stat=alloc_stat)
             if(alloc_stat.ne.0) call error_handler &
                  ("SEND_EIGVEC_ROT: deallocation 1 failed")
          enddo
          deallocate(eigvec_rot_real,eigvec_rot_imag,stat=alloc_stat)
          if(alloc_stat.ne.0) call error_handler &
               ("SEND_EIGVEC_ROT: deallocation 2 failed")
       endif
    else ! options_spin_orbit

       ! send n_rot to each slave
       if (send) then
          call comm_init_send(comm_all_other_hosts,msgtag_rot_levels)
          call commpack(n_rot(1,1),ssym%n_irrep*ssym%n_spin,1,info)
          if(info.ne.0) call error_handler &
               ("SEND_EIGVEC_ROT: packing failed")
       endif

       ! now copy those eigenvectors into eigvec_occ and eigvec_rot
       ! which have been changed by the perturbation theory and send
       ! eigvec_rot to each slave
       do i=1,ssym%n_irrep
          if (send) then
             allocate(eigvec_rot(ssym%dim(i),2,ssym%n_spin),stat=alloc_stat)
             if(alloc_stat.ne.0) call error_handler &
                  ("SEND_EIGVEC_ROT: allocation failed")
          endif
          do is=1,ssym%n_spin
             n=n_rot(i,is)
             if (n>0.and.n<=n_occo(is,i)) then
                eigvec_occ(i)%m(:,n,is) = eigvec(i)%m(:,n,is)
                if (n<n_occo(is,i)) then
                   eigvec_occ(i)%m(:,n+1,is) = eigvec(i)%m(:,n+1,is)
                endif
                if (send) then
                   eigvec_rot(:,1,is) = eigvec(i)%m(:,n,is)
                   if (n<n_occo(is,i)) then
                      eigvec_rot(:,2,is) = eigvec(i)%m(:,n+1,is)
                   else
                      eigvec_rot(:,2,is) = 0.0_r8_kind
                   endif
                endif
             endif
          enddo
          if (send) then
             length=ssym%dim(i)*2*ssym%n_spin
             if (length > 0) then
                call commpack(eigvec_rot(1,1,1),length,1,info)
                if(info.ne.0) call error_handler &
                     ("SEND_EIGVEC_ROT: packing failed")
             endif
             deallocate(eigvec_rot,stat=alloc_stat)
             if(alloc_stat.ne.0) call error_handler &
                  ("SEND_EIGVEC_ROT: deallocation failed")
          endif
       enddo
       if (send) then
          call comm_send()
       endif
    endif! options_spin_orbit
  end subroutine send_eigvec_rot

  !*************************************************************

  subroutine receive_eigvec_rot
    use comm_module
    use options_module, only: options_spin_orbit
    use occupation_module, only: n_occo, n_rot
    implicit none
    !** End of interface *****************************************

    integer(kind=i4_kind) :: i,is,n,ispin,alloc_stat,info,length
    real(kind=r8_kind), allocatable  :: eigvec_rot(:,:,:)
    type(arrmat2),allocatable         :: eigvec_rot_real(:),eigvec_rot_imag(:)

    integer(kind=i4_kind)                :: n_irrep
    integer(kind=i4_kind),allocatable    :: dim_irrep(:)
    ! n_irrep    : number of irreps
    ! dim_irrep  : number of independent functions in irrep

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
    ispin=ssym%n_spin

    ! get n_rot from the master
    allocate(n_rot(n_irrep,ispin),stat=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         ('receive_eigvec_rot: allocation (1) failed')
    call communpack(n_rot(1,1),n_irrep*ispin,1,info)
    if(info/=0) then
       call error_handler('receive_eigvec_rot: unpacking failed')
    endif

    if (options_spin_orbit) then
       !
       ! SPIN ORBIT
       !
       allocate(eigvec_rot_real(n_irrep),eigvec_rot_imag(n_irrep),stat=alloc_stat)
       if (alloc_stat/=0) call error_handler&
            ('receive_eigvec_rot: allocation (1) failed')
       do i=1,n_irrep
          allocate(eigvec_rot_real(i)%m(dim_irrep(i),2),eigvec_rot_imag(i)%m(dim_irrep(i),2),stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_rot: allocation (3) failed')
       enddo

       do i=1,n_irrep
          length=dim_irrep(i)*2
          if(length > 0) then
             call communpack(eigvec_rot_real(i)%m(1,1),length,1,info)
             call communpack(eigvec_rot_imag(i)%m(1,1),length,1,info)
          end if
       enddo

       ! Now update eigvec_occ and deallocate eigvec_rot again
       do i=1,n_irrep
          n=n_rot(i,1)
          if (n>0 .and. n<=n_occo(1,i)) then
             eigvec_occ_real(i)%m(:,n) = eigvec_rot_real(i)%m(:,1)
             eigvec_occ_imag(i)%m(:,n) = eigvec_rot_imag(i)%m(:,1)
             if (n<n_occo(1,i)) then
                eigvec_occ_real(i)%m(:,n+1) = eigvec_rot_real(i)%m(:,2)
                eigvec_occ_imag(i)%m(:,n+1) = eigvec_rot_imag(i)%m(:,2)
             endif
          endif
          deallocate( eigvec_rot_real(i)%m,eigvec_rot_imag(i)%m )
       enddo
       deallocate( n_rot, eigvec_rot_real,eigvec_rot_imag )
    else ! options_spin_orbit
       ! Now update eigvec_occ and deallocate eigvec_rot again
       do i=1,n_irrep
          allocate(eigvec_rot(ssym%dim(i),2,ispin),stat=alloc_stat)
          if (alloc_stat/=0) call error_handler&
               ('receive_eigvec_rot: allocation (2) failed')
          length=ssym%dim(i)*2*ispin
          if(length > 0) then
             call communpack(eigvec_rot(1,1,1),length,1,info)
          end if
          do is=1,ispin
             n=n_rot(i,is)
             if (n>0 .and. n<=n_occo(is,i)) then
                eigvec_occ(i)%m(:,n,is) = eigvec_rot(:,1,is)
                if (n<n_occo(is,i)) then
                   eigvec_occ(i)%m(:,n+1,is) = eigvec_rot(:,2,is)
                endif
             endif
          enddo
          deallocate( eigvec_rot )
       enddo
       deallocate( n_rot )
    endif! options_spin_orbit
    deallocate(dim_irrep)
    end subroutine receive_eigvec_rot

    ! ***************************************************************
!--------------- End of module ----------------------------------
 end module occupied_levels_module

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
module symm_adapt_module
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
       & RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !=======================================================
  !
  ! "Transition Density Matrix" struct:
  !

  type,private :: symm_adapt_pckdmat
     !
     ! packed type(symm_adapt_dmat), see below first
     !
     integer(IK)         :: n
     integer(IK),pointer :: m2(:),m1(:) ! m(n)
     real(RK),pointer    :: c(:)        ! c(n)

     !-------- private ----------------
     integer(IK)               :: L2,L1
     integer(IK),pointer       :: ibuf(:,:)
     !
     ! ibuf((2*L2+1) * (2*L1+1), 2)
     !
     real(RK),pointer          :: rbuf(:)
     !
     ! rbuf((2*L2+1) * (2*L1+1))
     !
     ! m2 => ibuf(1:n,2)
     ! m1 => ibuf(1:n,1)
     ! c  => rbuf(1:n)
  end type symm_adapt_pckdmat
  
  type,private :: symm_adapt_dmat
     !
     ! Rho(m2,m1) = SUM_a,p { C*2(a,p,m2) * C1(a,p,m1) }
     ! (a,p) are spinor component and partner of an irrep
     ! such a construct is used in symmetry adaption procedure

     integer(IK)               :: L2,L1

     real(RK),pointer          :: dmat(:,:,:)
     !
     ! dmat(2*L2+1,2*L1+1,2) - real and imaginary parts
     logical,pointer           :: mask(:,:,:)
     type(symm_adapt_pckdmat)  :: pck(2) 
     ! re and im
  end type symm_adapt_dmat


  !------------ Interface statements ------------------------------
 
  interface alloc
     module procedure alloc_symm_adapt_dmat
     module procedure alloc_symm_adapt_pckdmat
  end interface

  interface dealloc
     module procedure dealloc_symm_adapt_dmat
     module procedure dealloc_symm_adapt_pckdmat
  end interface

  interface free      ! recursive dealloc
     module procedure dealloc_symm_adapt_dmat
  end interface

  interface tpack
     module procedure pack_symm_adapt_dmat
  end interface

  interface CalcDim
     module procedure CalcDim_one
     module procedure CalcDim_two
  end interface

  interface tload
     module procedure load_ua_to_symm_adapt_ul
  end interface

  !------------ public functions and subroutines ------------------

  public :: symm_adapt_module_setup!()
  public :: symm_adapt_module_close!()

  public :: symm_adapt_2c
  public :: symm_adapt_3c
  public :: symm_adapt_2cv
  public :: symm_adapt_3cv

  !------------ Declaration of constants and variables ----


  logical,parameter :: YES = .true., NO = .false.

  logical,private   :: debug       = YES
  logical,private   :: initialized = NO

  logical,private   :: SpinOrbit
  logical,private   :: Small

  !------------ Subroutines ---------------------------------------
contains

  subroutine symm_adapt_module_setup()
    use error_module
    use comm_module
    use options_module,    only: options_spin_orbit
    use symm_adapt_struct, only: init_symm_adapt_struct
    use symm_adapt_xpack
    use spin_orbit_module, only: is_on, op_FitTrafo
    implicit none
    ! *** end of interface ***

    DPRINT MyID, 'sam/symm_adapt_module_setup: entered'

    if ( initialized ) then
        WARN('why calling?')
        return
    endif

    SpinOrbit = options_spin_orbit

    Small = is_on(op_FitTrafo)

    if(.not.comm_parallel())then
       !
       ! serial run:
       !
       DPRINT MyID,'sam/symm_adapt_module_setup: serial run'

       call PreAlloc()
       call SymAdpLoad()
       call Complete()
    else
       !
       ! parallel run:
       !
       DPRINT MyID,'sam/symm_adapt_module_setup: parallel run'

       call PreAlloc()

       if(comm_i_am_master())then
          !
          ! master:
          !
          call SymAdpLoad()

          DPRINT MyID,'sam/symm_adapt_module_setup: trying to send ...'
          call symm_adapt_send()
       else
          !
          ! slaves:
          !
          DPRINT MyID,'sam/symm_adapt_module_setup: trying to receive ...'
          call symm_adapt_receive()
       endif

       call Complete()
    endif

    call init_symm_adapt_struct()

    initialized = YES
    DPRINT MyID,'sam/symm_adapt_module_setup: exited'
  contains

    subroutine PreAlloc()
      use unique_atom_module, only: n_ua => n_unique_atoms
      use symm_adapt_struct
      use error_module
      implicit none
      ! *** end of interface ***

      integer :: memstat

      !---------------------------------
      !
      ! do allocation >>>
      !
      DPRINT MyID,'sam/symm_adapt_module_setup: doing InitPreAlloc'

      if(.not.SpinOrbit)then
!        WARN('PreAlloc: no SO, ignore')
         return
      endif

      allocate(&
           & LSymAdp(n_ua), &
           & LSymAdpL(n_ua),&
           & LSymAdpS(n_ua),&
           & STAT=memstat)
      call error(memstat,"sam/symm_adapt_module_setup: alloc 1 failed")
    end subroutine PreAlloc

    subroutine SymAdpLoad()
      use error_module
      use uatom_symmadapt, only: uaSymmSpor, uas_Large, uas_Small
      use symm_adapt_struct
      implicit none
      ! *** end of interface ***

      integer(IK) :: a,n_ua

      !---------------------------------
      !
      ! load symm-adapt from ua-module >>>
      !
      DPRINT MyID,'sam/symm_adapt_module_setup: master doing loading...'

      if(.not.SpinOrbit)then
!        WARN('SymAdpLoad: no SO, ignore')
         return
      endif

      n_ua = size(uaSymmSpor)

      DPRINT MyID,'sam/symm_adapt_module_setup: loading uaSymm ...'
      do a=1,n_ua
         call tload(uaSymmSpor(a),LSymAdp(a))
      enddo

      DPRINT MyID,'sam/symm_adapt_module_setup: loading uas_Large ...'
      do a=1,n_ua
!!$         DPRINT 'sam/symm_adapt_module_setup: ua=',a
         call tload(uas_Large(a),LSymAdpL(a))
!!$         call show_ul(LSymAdpL(a))
      enddo

      if(Small)then
         DPRINT MyID,'sam/symm_adapt_module_setup: loading uas_Small ...'
         do a=1,n_ua
!!$            DPRINT 'sam/symm_adapt_module_setup: ua=',a
            call tload(uas_Small(a),LSymAdpS(a))
!!$            call show_ul(LSymAdpS(a))
         enddo
      endif
    end subroutine SymAdpLoad

    subroutine Complete()
      use error_module
      use unique_atom_module, only: n_ua => n_unique_atoms
      use dimensions
      use symm_adapt_struct
      use uatom_symmadapt, only: uaSymm
      implicit none
      ! *** end of interface ***

      integer(IK) :: a

      !---------------------------------
      !
      ! complete initialization >>>
      !
      DPRINT MyID,'sam/symm_adapt_module_setup: doing Complete'
      
      DPRINT MyID,'sam/symm_adapt_module_setup: call dimensions_init(InitLevelAlloc)'
      call dimensions_init(InitLevelAlloc)

      DPRINT MyID,'sam/symm_adapt_module_setup: Load vector dimensions ...'
      ASSERT(n_ua==size(uaSymm))
      ASSERT(n_ua==size(SymDim))
      do a=1,n_ua
         call CalcDim(uaSymm(a), SymDim(a))
      enddo

      if(.not.SpinOrbit)then
!        WARN('Complete: no SO, skip')
         goto 999
      endif

      DPRINT MyID,'sam/symm_adapt_module_setup: Load projective dimensions ...'
      do a=1,n_ua
         call CalcDim(LSymAdp(a), SymDimSpor(a))
      enddo

      DPRINT MyID,'sam/symm_adapt_module_setup: Load projective dimensions (Large)...'
      do a=1,n_ua
         call CalcDim(LSymAdpL(a), SymDimSporL(a))
      enddo

      if(Small)then
         DPRINT MyID,'sam/symm_adapt_module_setup: Load projective dimensions (Small)...'
         do a=1,n_ua
            call CalcDim(LSymAdpS(a), SymDimSporS(a))
         enddo
      endif

999   continue ! jump here
      DPRINT MyID,'sam/symm_adapt_module_setup: call dimensions_init(InitLevelComplete)'
      call dimensions_init(InitLevelComplete)
      
    end subroutine Complete

  end subroutine symm_adapt_module_setup

  subroutine symm_adapt_module_close()
    use error_module
    use symm_adapt_struct
    use dimensions, only: dimensions_free
    use spin_orbit_module, only: is_on, op_FitTrafo
    implicit none
    ! *** end of interface ***

    integer :: memstat

    DPRINT MyID, "sas/symm_adapt_module_close: entered"

    if ( .not.initialized ) then
        WARN('why calling?')
        return
    endif

    if(.not.SpinOrbit)then
!      WARN('symm_adapt_module_close: no SO, skip')
       goto 999
    endif

    call free(LSymAdp)
    call free(LSymAdpL)
    if(is_on(op_FitTrafo))then
       call free(LSymAdpS)
    endif

    deallocate(&
         & LSymAdp,LSymAdpL,LSymAdpS,&
         & STAT=memstat)
    call error(memstat,"sas/symm_adapt_module_close: dealloc failed")

999 continue ! jump here
    DPRINT MyID,"sas/symm_adapt_module_close: call dimensions_free()"
    call dimensions_free()

    initialized = NO
  end subroutine symm_adapt_module_close

  !=======================================================

  subroutine CalcDim_one(sa, d)
    use dimensions, only: SubSpaceDim, alloc
    use symm_adapt_struct
    implicit none
    type(symm_adapt_ul),intent(in)  :: sa
    type(SubSpaceDim),intent(inout) :: d
    ! *** end of interface ***

    integer(IK) :: n_irr,lmax,l,irr

    n_irr = sa%n_irr
    lmax  = sa%lmax

    call alloc(lmax,n_irr,d)

    l_: do l=0,lmax
       d%LM(l,:) = 0
       irr_: do irr=1,n_irr
          if(sa%exists(l,irr))then
             d%LM(l,irr) = d%LM(l,irr) + sa%aos(1,l,irr)%n_indep
          endif
       enddo irr_
    enddo l_
  end subroutine CalcDim_one

  subroutine CalcDim_two(ua, d)
    use error_module
    use dimensions, only: SubSpaceDim, alloc
    use symm_adapt_struct
    use uatom_symmadapt
    implicit none
    type(uatom),intent(in)          :: ua
    type(SubSpaceDim),intent(inout) :: d
    ! *** end of interface ***

    integer(IK) :: l,irr,lmax,n_irr

    if(ua%spor)call error("sam/CalcDim_two: use CalcDim_one")

    n_irr = ua%n_irr
    lmax  = ua%lmax

    call alloc(lmax,n_irr,d)

    l_: do l=0,lmax
       d%LM(l,:) = 0
       irr_: do irr=1,n_irr
          d%LM(l,irr) = d%LM(l,irr)&
               & + ua%symadapt_partner(irr,l)%n_independent_fcts
       enddo irr_
    enddo l_
  end subroutine CalcDim_two

  subroutine load_ua_to_symm_adapt_ul(ua,sa)
    !
    ! allocates memory inbetween
    !
    use uatom_symmadapt, only: uatom, partner_type, sa_int_type
    use error_module
    use symm_adapt_struct
    implicit none
    type(uatom),intent(in)            :: ua
    type(symm_adapt_ul),intent(inout) :: sa
    ! ***  end of interface ***

    integer(IK) :: n_irr,irr,lmax,l,n_ea,e,n_indep,f,dim_irr,p,c

    type(symm_adapt_irrbas),pointer        :: irrbas
    type(symm_adapt_ao),pointer            :: o

    type(sa_int_type) :: si

    n_irr =   size(ua%symadapt_spor_partner,1)
    lmax  = ubound(ua%symadapt_spor_partner,2)

    DPRINT 'load_ua_to_symm_adapt_ua: lmax=',lmax,' n_irr=',n_irr

    ASSERT(ua%lmax==lmax)
    ASSERT(ua%n_irr==n_irr)

    n_ea  =   n_equiv_atoms(ua)

    call alloc(n_ea,lmax,n_irr,sa)

    l_: do l=0,lmax
       irr_: do irr=1,n_irr

          n_indep = ua%symadapt_spor_partner(irr, l)%n_independent_fcts

          if(n_indep>0)then

             sa%exists(l,irr) = .true.

             e_: do e = 1, n_ea

                dim_irr = size(ua%symadapt_spor_partner(irr, l)%sa_spor_int, 4)

                irrbas => sa%aos(e, l, irr)

                call alloc(dim_irr, n_indep, irrbas)

                p_: do p = 1, dim_irr
                   f_: do f = 1, n_indep
                      c_: do c = 1, 2

                         o  => irrbas%bas(p,f)%psi(c)
                         ! FIXME: a copy or a ref?
                         si = ua%symadapt_spor_partner(irr, l)%sa_spor_int(e, c, f, p)

                         call alloc(si%n_fcts, o)

                         o%m  = si%m
                         o%re = si%re
                         o%im = si%im
                      enddo c_
                   enddo f_
                enddo p_
             enddo e_
          else
             sa%exists(l,irr) = .false.
          endif
       enddo irr_
    enddo l_
  end subroutine load_ua_to_symm_adapt_ul

  !********************************************
  !
  !      SYMMETRY ADAPTION PRECEDURES:
  !
  !********************************************

  subroutine symm_adapt_2c( p2c, sa2c, &
                            ea2,L2,sa2,&
                            ea1,L1,sa1,&
                            weight )
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use symm_adapt_struct, only: symm_adapt_ul, symm_adapt_irrbas
    use symm_adapt_int, only: sa_int_block
    implicit none
    external integral_interrupt_2cob3c
    integer(IK),         intent(in)     :: L1,L2,ea1,ea2
    real(RK),            intent(in)     :: p2c(:,:,:,:)
    type(sa_int_block),  intent(inout)  :: sa2c(:) ! (n_irr)
    type(symm_adapt_ul), intent(in)     :: sa2,sa1 ! (n_irr)
    real(RK),            intent(in)     :: weight
    ! *** end of interface ***

    type(sa_int_block)              :: b2c

    type(symm_adapt_irrbas),pointer :: irrb1,irrb2

    integer(IK) :: irr
    integer(IK) :: f1,f2
    real(RK)    :: weight2

    type(symm_adapt_dmat) :: dmat

    integer(IK)                      :: i,n
    integer(IK),pointer,dimension(:) ::  m2, m1
    real(RK),pointer,dimension(:)    :: cs
   
    call alloc(L2,L1,dmat)

    ASSERT(sa1%n_irr==sa2%n_irr)

    irr_: do irr=1,sa1%n_irr

       if(.not.(sa1%exists(L1,irr).and.sa2%exists(L2,irr))) cycle

       irrb1 => sa1%aos(ea1,L1,irr)
       irrb2 => sa2%aos(ea2,L2,irr)

       ASSERT(irrb1%dim_irr==irrb2%dim_irr)

       weight2 = weight/irrb1%dim_irr

       b2c  = sa2c(irr)

       f1_: do f1=1,irrb1%n_indep
          f2_: do f2=1,irrb2%n_indep


             call tdensmat(irrb2%bas(:,f2),irrb1%bas(:,f1),dmat)

             !---------------------------------
             ! two center integrals:
             !

             n   =  dmat%pck(1)%n
             cs  => dmat%pck(1)%c
              m2 => dmat%pck(1)%m2
              m1 => dmat%pck(1)%m1

             cs = cs * weight2

             do i=1,n
                b2c%re(:,:,f2,f1) = b2c%re(:,:,f2,f1)&
                     & +cs(i) * p2c(:,:, m2(i), m1(i))
             enddo

             n   =  dmat%pck(2)%n
             cs  => dmat%pck(2)%c
              m2 => dmat%pck(2)%m2
              m1 => dmat%pck(2)%m1

             cs = cs * weight2

             do i=1,n
                b2c%im(:,:,f2,f1) = b2c%im(:,:,f2,f1)&
                     & +cs(i) * p2c(:,:, m2(i), m1(i))
             enddo


          enddo f2_
       enddo f1_

    enddo irr_

    call free(dmat)
  end subroutine symm_adapt_2c

  subroutine symm_adapt_3c( p3c, sa3c,&
                            ea2,L2,sa2,&
                            ea1,L1,sa1,&
                            weight &
                            )
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use symm_adapt_struct, only: symm_adapt_ul, symm_adapt_irrbas
    use symm_adapt_int, only: sa_3c_int_block
    implicit none
    integer(IK),         intent(in)     :: L1,L2,ea1,ea2
    real(RK),            intent(in)     :: p3c(:,:,:,:,:)
    type(sa_3c_int_block), intent(inout) :: sa3c(:) ! (n_irr)
    type(symm_adapt_ul), intent(in)     :: sa2,sa1 ! (n_irr)
    real(RK),            intent(in)     :: weight
    ! *** end of interface ***

    real(RK),parameter              :: zero = 0.0_rk, one = 1.0_rk

    type(sa_3c_int_block)           :: b3c

    type(symm_adapt_irrbas),pointer :: irrb1,irrb2

    integer(IK) :: irr
    integer(IK) :: f1,f2
    real(RK)    :: weight2

    type(symm_adapt_dmat) :: dmat

    integer(IK)                      :: i,n
    integer(IK),pointer,dimension(:) ::  m2, m1
    real(RK),pointer,dimension(:)    :: cs
   
    call alloc(L2,L1,dmat)

    ASSERT(sa1%n_irr==sa2%n_irr)

    irr_: do irr=1,sa1%n_irr

       if(.not.(sa1%exists(L1,irr).and.sa2%exists(L2,irr))) cycle

       irrb1 => sa1%aos(ea1,L1,irr)
       irrb2 => sa2%aos(ea2,L2,irr)

       ASSERT(irrb1%dim_irr==irrb2%dim_irr)

       weight2 = weight/irrb1%dim_irr

       b3c  = sa3c(irr)

       f1_: do f1=1,irrb1%n_indep
          f2_: do f2=1,irrb2%n_indep


             call tdensmat(irrb2%bas(:,f2),irrb1%bas(:,f1),dmat)

             !---------------------------------
             ! three center integrals:
             !
             n   =  dmat%pck(1)%n
             cs  => dmat%pck(1)%c
              m2 => dmat%pck(1)%m2
              m1 => dmat%pck(1)%m1

             cs = cs * weight2

             do i=1,n
                b3c%re(:,:,:,f2,f1) = b3c%re(:,:,:,f2,f1)&
                     & + cs(i) * p3c(:,:,:, m2(i), m1(i))
             enddo

             n   =  dmat%pck(2)%n
             cs  => dmat%pck(2)%c
              m2 => dmat%pck(2)%m2
              m1 => dmat%pck(2)%m1

             cs = cs * weight2

             do i=1,n 
                b3c%im(:,:,:,f2,f1) = b3c%im(:,:,:,f2,f1)&
                     & + cs(i) * p3c(:,:,:, m2(i), m1(i))
             enddo

          enddo f2_
       enddo f1_

    enddo irr_

    call free(dmat)
  end subroutine symm_adapt_3c

  subroutine symm_adapt_2cv( p2cv,sa2cv,&
                             ea2,L2,sa2,&
                             ea1,L1,sa1,&
                             weight )
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use symm_adapt_struct, only: symm_adapt_ul, symm_adapt_irrbas
    use symm_adapt_int, only: sa_int_block
    implicit none
    external integral_interrupt_2cob3c
    integer(IK),         intent(in)     :: L1,L2,ea1,ea2
    real(RK),            intent(in)     :: p2cv(:,:,:,:,:)
    type(sa_int_block),  intent(inout)  :: sa2cv(:) ! (n_irr)
    type(symm_adapt_ul), intent(in)     :: sa2,sa1 ! (n_irr)
    real(RK),            intent(in)     :: weight
    ! *** end of interface ***

    real(RK),parameter              :: zero = 0.0_rk, one = 1.0_rk

    type(sa_int_block)              :: b2cv

    type(symm_adapt_irrbas),pointer :: irrb1,irrb2

    integer(IK) :: irr
    integer(IK) :: f1,f2,xyz
    real(RK)    :: weight2

    type(symm_adapt_dmat) :: sdmat(3)

    integer(IK)                      :: i,n
    integer(IK),pointer,dimension(:) ::  m2, m1
    real(RK),pointer,dimension(:)    :: cs
   
    do i=1,3
       call alloc(L2,L1,sdmat(i))
    enddo

    ASSERT(sa1%n_irr==sa2%n_irr)

    irr_: do irr=1,sa1%n_irr

       if(.not.(sa1%exists(L1,irr).and.sa2%exists(L2,irr))) cycle

       irrb1 => sa1%aos(ea1,L1,irr)
       irrb2 => sa2%aos(ea2,L2,irr)

       ASSERT(irrb1%dim_irr==irrb2%dim_irr)

       weight2 = weight/irrb1%dim_irr

       b2cv = sa2cv(irr)

       f1_: do f1=1,irrb1%n_indep
          f2_: do f2=1,irrb2%n_indep


             call tdensmat(irrb2%bas(:,f2),irrb1%bas(:,f1),S=sdmat)

             !---------------------------------
             !
             ! do sigmas
             !
             xyz_: do xyz=1,3

                n   =  sdmat(xyz)%pck(1)%n
                cs  => sdmat(xyz)%pck(1)%c
                 m2 => sdmat(xyz)%pck(1)%m2
                 m1 => sdmat(xyz)%pck(1)%m1

                cs = cs * weight2

                do i=1,n
                   b2cv%re(:,:,f2,f1) = b2cv%re(:,:,f2,f1)&
                        & +cs(i) * p2cv(:,:, m2(i), m1(i),xyz)
                enddo

                n   =  sdmat(xyz)%pck(2)%n
                cs  => sdmat(xyz)%pck(2)%c
                 m2 => sdmat(xyz)%pck(2)%m2
                 m1 => sdmat(xyz)%pck(2)%m1

                cs = cs * weight2

                do i=1,n
                   b2cv%im(:,:,f2,f1) = b2cv%im(:,:,f2,f1)&
                        & +cs(i) * p2cv(:,:, m2(i), m1(i),xyz)
                enddo
             enddo xyz_
             !
             ! done sigmas
             !---------------------------------

          enddo f2_
       enddo f1_

    enddo irr_

    do i=1,3
       call free(sdmat(i))
    enddo
  end subroutine symm_adapt_2cv

  subroutine symm_adapt_3cv( p3cv, sa3cv,&
                             ea2,L2,sa2,&
                             ea1,L1,sa1,&
                             weight )
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use symm_adapt_struct, only: symm_adapt_ul, symm_adapt_irrbas
    use symm_adapt_int, only: sa_3c_int_block
    implicit none
    integer(IK),         intent(in)     :: L1,L2,ea1,ea2
    real(RK),            intent(in)     :: p3cv(:,:,:,:,:,:)
    type(sa_3c_int_block), intent(inout) :: sa3cv(:) ! (n_irr)
    type(symm_adapt_ul), intent(in)     :: sa2,sa1 ! (n_irr)
    real(RK),            intent(in)     :: weight
    ! *** end of interface ***

    real(RK),parameter              :: zero = 0.0_rk, one = 1.0_rk

    type(sa_3c_int_block)           :: b3cv

    type(symm_adapt_irrbas),pointer :: irrb1,irrb2

    integer(IK) :: irr
    integer(IK) :: f1,f2,xyz
    real(RK)    :: weight2

    type(symm_adapt_dmat) :: sdmat(3)

    integer(IK)                      :: i,n
    integer(IK),pointer,dimension(:) ::  m2, m1
    real(RK),pointer,dimension(:)    :: cs
   
    do i=1,3
       call alloc(L2,L1,sdmat(i))
    enddo

    ASSERT(sa1%n_irr==sa2%n_irr)

    irr_: do irr=1,sa1%n_irr

       if(.not.(sa1%exists(L1,irr).and.sa2%exists(L2,irr))) cycle

       irrb1 => sa1%aos(ea1,L1,irr)
       irrb2 => sa2%aos(ea2,L2,irr)

       ASSERT(irrb1%dim_irr==irrb2%dim_irr)

       weight2 = weight/irrb1%dim_irr

       b3cv  = sa3cv(irr)

       f1_: do f1=1,irrb1%n_indep
          f2_: do f2=1,irrb2%n_indep


             call tdensmat(irrb2%bas(:,f2),irrb1%bas(:,f1),S=sdmat)

             !---------------------------------
             !
             ! do sigmas
             !
             xyz_: do xyz=1,3

                n   =  sdmat(xyz)%pck(1)%n
                cs  => sdmat(xyz)%pck(1)%c
                 m2 => sdmat(xyz)%pck(1)%m2
                 m1 => sdmat(xyz)%pck(1)%m1

                cs = cs * weight2

                do i=1,n
                   b3cv%re(:,:,:,f2,f1) = b3cv%re(:,:,:,f2,f1)&
                        & +cs(i) * p3cv(:,:,:, m2(i), m1(i),xyz)
                enddo

                n   =  sdmat(xyz)%pck(2)%n
                cs  => sdmat(xyz)%pck(2)%c
                 m2 => sdmat(xyz)%pck(2)%m2
                 m1 => sdmat(xyz)%pck(2)%m1

                cs = cs * weight2

                do i=1,n
                   b3cv%im(:,:,:,f2,f1) = b3cv%im(:,:,:,f2,f1)&
                        & +cs(i) * p3cv(:,:,:, m2(i), m1(i),xyz)
                enddo
             enddo xyz_
             !
             ! done sigmas
             !---------------------------------

          enddo f2_
       enddo f1_

    enddo irr_

    do i=1,3
       call free(sdmat(i))
    enddo
  end subroutine symm_adapt_3cv

  subroutine tdensmat(b2,b1,d,s) 
    !
    ! Transition Density Matrix, dont take the name seriously
    !
    use symm_adapt_struct
    implicit none
    type(symm_adapt_spinor),intent(in)  :: b2(:),b1(:) ! b(dim_irr)
    type(symm_adapt_dmat),intent(inout) :: d,s(3)
    optional                            :: d,s
    target :: b1,b2,d
    ! *** end of interface ***

    real(RK),parameter               :: zero = 0.0_rk
    integer(IK)                      :: p,c,i1,i2,m1,m2,r,k,xyz
    integer(IK)                      :: sa,sb,a,b,cf
    type(symm_adapt_ao),pointer      :: o1,o2
    real(RK),pointer                 :: dmat(:,:,:)

    if(present(d))then
       dmat => d%dmat
       dmat = zero

       do p=1,size(b1)
          do c=1,2

             o1 => b1(p)%psi(c)
             o2 => b2(p)%psi(c)

             do i1=1,o1%n
                m1 = o1%m(i1)

                do i2=1,o2%n
                   m2 = o2%m(i2)

                   dmat(m2,m1,1) = dmat(m2,m1,1)&
                        & + o2%re(i2) * o1%re(i1) + o2%im(i2) * o1%im(i1)

                   dmat(m2,m1,2) = dmat(m2,m1,2)&
                        & - o2%im(i2) * o1%re(i1) + o2%re(i2) * o1%im(i1)

                enddo
             enddo
          enddo
       enddo

       call tpack(d)
    endif

    if(present(s))then
       xyz_:do xyz=1,3

          dmat => s(xyz)%dmat
          dmat =  zero

          r_:do r=1,2
             p_:do p=1,size(b1)
                k_:do k=1,4

                   sa = tSIGMA(r,xyz)%alpha(k)
                   sb = tSIGMA(r,xyz)%beta(k)

                   a  = tSIGMA(r,xyz)%a(k)
                   b  = tSIGMA(r,xyz)%b(k)
                   cf = tSIGMA(r,xyz)%c(k)

                   o1 => b1(p)%psi(sb)
                   o2 => b2(p)%psi(sa)

                   i1_:do i1=1,o1%n
                      m1 = o1%m(i1)
                      i2_:do i2=1,o2%n
                         m2 = o2%m(i2)

                         dmat(m2,m1,r) = dmat(m2,m1,r)&
                              & + cf * o2%c(i2,a) * o1%c(i1,b)

                      enddo i2_
                   enddo i1_
                enddo k_
             enddo p_
          enddo r_

          ! i-factor >>>
          dmat(:,:,1:2) =   dmat(:,:,2:1:-1)
          dmat(:,:,  1) = - dmat(:,:,1) 
       enddo xyz_

       !debug>>>
       s(1)%dmat = - s(1)%dmat  !<<<HERE!!!
       s(2)%dmat = - s(2)%dmat  !<<<HERE!!!
       s(3)%dmat = + s(3)%dmat
       !<<<debug

       do xyz=1,3
          call tpack(s(xyz))
       enddo
    endif
  end subroutine tdensmat

  subroutine pack_symm_adapt_dmat(d)
    implicit none
    type(symm_adapt_dmat),intent(inout) :: d
    target :: d
    ! *** end of interface ***

    real(RK),parameter :: zero = 0.0_rk
    integer(IK)        :: n,r,i
    type(symm_adapt_pckdmat),pointer :: dp
    logical,pointer  :: mask(:,:,:)
    real(RK)         :: eps

    eps = epsilon(zero) * 10000.0_rk

    mask => d%mask

    do r=1,2 ! re an im
       dp => d%pck(r)

       mask(:,:,r) = (abs(d%dmat(:,:,r)) > eps)
       ! (abs(d%dmat(:,:,r)) /= zero ) 
       !.true. ! (abs(dmat(:,:,r)) > eps)
       
       n     =  count(mask(:,:,r))
       dp%n  =  n
       dp%c  => dp%rbuf(1:n)
       dp%m2 => dp%ibuf(1:n,2)
       dp%m1 => dp%ibuf(1:n,1)

       if(n==0)cycle !<<< for HP

       dp%c  = pack(d%dmat(:,:,r),mask(:,:,r))
       dp%m2 = pack(spread((/(i,i=1,2*dp%L2+1)/), 2, 2*dp%L1+1),mask(:,:,r))
       dp%m1 = pack(spread((/(i,i=1,2*dp%L1+1)/), 1, 2*dp%L2+1),mask(:,:,r))
    enddo
  end subroutine pack_symm_adapt_dmat

  !-------------------------------------------------------
  !
  ! MEMORY MANAGMENT PROCEDURES
  !

  subroutine alloc_symm_adapt_dmat(L2,L1,d)
    use error_module
    implicit none
    integer(IK),intent(in)              :: L2,L1
    type(symm_adapt_dmat),intent(inout) :: d
    ! *** end of interface ***

    integer     :: memstat
    integer(IK) :: n_lm1,n_lm2

    d%L1 = L1
    d%L2 = L2

    n_lm1 = 2*L1+1
    n_lm2 = 2*L2+1

    allocate(d%dmat(n_lm2,n_lm1,2),d%mask(n_lm2,n_lm1,2),STAT=memstat)
    call error(memstat,"sam/alloc_symm_adapt_dmat: alloc failed")

    call alloc(L2,L1,d%pck(1))
    call alloc(L2,L1,d%pck(2))
  end subroutine alloc_symm_adapt_dmat

  subroutine dealloc_symm_adapt_dmat(d)
    use error_module
    implicit none
    type(symm_adapt_dmat),intent(inout) :: d
    ! *** end of interface ***

    integer     :: memstat

    d%L1 = -1
    d%L2 = -1

    deallocate(d%dmat,d%mask,STAT=memstat)
    call error(memstat,"sam/dealloc_symm_adapt_dmat: alloc failed")

    call dealloc(d%pck(1))
    call dealloc(d%pck(2))
  end subroutine dealloc_symm_adapt_dmat

  subroutine alloc_symm_adapt_pckdmat(L2,L1,d)
    use error_module
    implicit none
    integer(IK),intent(in)                 :: L2,L1
    type(symm_adapt_pckdmat),intent(inout) :: d
    ! *** end of interface ***

    integer     :: memstat
    integer(IK) :: n_lm1,n_lm2

    d%L1 = L1
    d%L2 = L2

    n_lm1 = 2*L1+1
    n_lm2 = 2*L2+1

    allocate(d%ibuf(n_lm2*n_lm1,2), d%rbuf(n_lm2*n_lm1),STAT=memstat)
    call error(memstat,"sam/alloc_symm_adapt_pckdmat: alloc failed")

    d%n = 0
    d%m2 => d%ibuf(1:d%n,2)
    d%m1 => d%ibuf(1:d%n,1)
    d%c  => d%rbuf(1:d%n)
  end subroutine alloc_symm_adapt_pckdmat

  subroutine dealloc_symm_adapt_pckdmat(d)
    use error_module
    implicit none
    type(symm_adapt_pckdmat),intent(inout) :: d
    ! *** end of interface ***

    integer     :: memstat

    d%L1 = -1
    d%L2 = -1

    deallocate(d%ibuf, d%rbuf,STAT=memstat)
    call error(memstat,"sam/dealloc_symm_adapt_pckdmat: dealloc failed")

    d%n = -1
    nullify(d%m2,d%m1)
    nullify(d%c)
  end subroutine dealloc_symm_adapt_pckdmat

  !
  ! MISCELANEOUS:
  !

  function n_equiv_atoms(a) result(n)
    use uatom_symmadapt, only: uatom, partner_type
    implicit none
    type(uatom),intent(in) :: a
    integer(IK) :: n
    ! *** end of interface ***

    integer(IK) :: irr, L, n_ea
    logical     :: found

    found = .false.

    do irr=1,size(a%symadapt_spor_partner,1)
       do L=0,ubound(a%symadapt_spor_partner,2)

          if(a%symadapt_spor_partner(irr, L)%n_independent_fcts > 0 )then

             n_ea = size(a%symadapt_spor_partner(irr, L)%sa_spor_int, 1)

             if(found)then
                ASSERT(n==n_ea)
             else
                n = n_ea
                found = .true.
             endif
          endif
       enddo
    enddo
  end function n_equiv_atoms

  !--------------- End of module ----------------------------------
end module symm_adapt_module

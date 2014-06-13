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
module shgi_slv
  !---------------------------------------------------------------
  !
  ! Integrals due to the charges of the solvation shell.
  !
  ! Copyright (c) 2006-2008 Alexey Shor
  ! Copyright (c) 2007-2013 Alexei Matveev
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
! use CPU_TIME for timers:
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind, & ! type specification parameters
       I8K=>i8_kind
  use shgi_cntrl
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: shgi_gr_solv_drv
  public :: shgi_gr_solv_drv_vtn
  public :: shgi_gr_Q_solv_drv
  public :: shgi_sd_solv_drv

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  ! ALL INTEGER CONSTANTS, LIKE GAX,GAY,GAZ ARE NOW IN
  !                     shgi_cntrl.f90
  ! THIS HAS BEEN DONE TO ENABLE ITS USE IN OTHER MODULES

  ! MOST GLOBAL VARIABLES HOLDING ``ANGULAR'' FACTORS WERE MOVED TO
  !                     shgi_common.f90
  ! THIS HAS BEEN DONE TO SPLIT THIS FILE INTO PARTS LATER

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !****************************************************************
  !****************** DRIVER FOR GRADIENTS ************************
  !****************************************************************

  subroutine shgi_gr_solv_drv(IU1,IE1,IL1,IU2,IE2,IL2,uas,densmat,VECUAGR,VECPCGR,MATUAGR,MATPCGR)
    use unique_atom_module, only: uat=>unique_atom_type,n_unique_atoms
    use solv_cavity_module, only: to_calc_grads, with_pc, fixed_pc
    use pointcharge_module, only: pointcharge_N
    use datatype, only: arrmat4 ! PG storage of grs
    use options_module, only: options_integral_expmax
    use shgi_common, only: CUTOFF
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: shgi_gr_wd_store, shgi_timing, &
                           shgi_gr_abc_store, shgi_gr_wd_to_abc, &
                           shgi_ve_abc_store
    use shgi_pcm, only: shgi_gr_pc
    implicit none
    type(uat)  , intent(in)      :: uas(:)            ! array of unique atoms, normally all of them
    integer(IK), intent(in)      :: IU1,IE1,IL1,IU2,IE2,IL2
    real(RK)     , intent(in)    :: densmat(:,:,:,:)  ! (NEA,NEB,2*LA+1,2*LB+1), section of the density matrix
    real(RK)     , intent(inout) :: VECUAGR(:)        ! (num of tot-sym grads)            -- grads wrt UAs
    real(RK)     , intent(inout) :: VECPCGR(:)        ! (num of tot-sym grads)            -- grads wrt PCs
    type(arrmat4), intent(inout) :: MATUAGR(:)   ! (num of tot-sym grads)%m(:,:,:,:) -- matrix grads wrt UAs
    type(arrmat4), intent(inout) :: MATPCGR(:)   ! (num of tot-sym grads)%m(:,:,:,:) -- matrix grads wrt PCs
    optional :: VECUAGR, VECPCGR, densmat ! -- in "property" runs (forces) we dont need matrices
    optional :: MATUAGR, MATPCGR          ! -- in "pre-CPKS" run we need the matrices
    ! *** end of interface ***

    integer(IK) :: up,ep,N_points,n_equal,ua,N_pc
    real(RK)    :: x(3),z
    integer(IK) :: ism ! symmetry type???
    integer(IK) :: memstat
    integer(IK) :: ma,mb,i
    real(RK), allocatable :: GRSOLV(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6), W,D
    real(RK), allocatable :: GRSUM(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6), W,D
    real(RK), allocatable :: GAB(:,:,:,:)    ! (NAB,2*LA+1,2*LB+1,9) -- wrt A, B, and C
    real(RK), allocatable :: PAB(:,:,:)      ! (NAB,2*LA+1,2*LB+1)   -- packed densmat
    real(RK)              :: VAB(9)          !                   (9) -- wrt A, B, and C, trace(PAB,GAB)
    logical               :: have_densmat

    DPRINT  'SHGI: shgi_gr_solv_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

    FPP_TIMER_START(tot)
    FPP_TIMER_START(totG)
    FPP_TIMER_START(pcs)
    FPP_TIMER_START(slgr)

    !
    ! FIXME:  what is  wrong with  the value  of pointcharge_N  if the
    ! condition does not hold? Better make pointcharge_N always have a
    ! reasonable value,  zero or not,  irrespective of what  the flags
    ! are:
    !
    if (with_pc .and. .not. fixed_pc) then
       N_pc = pointcharge_N
    else
       N_pc = 0
    endif

    N_points = to_calc_grads%n_points

    have_densmat = present(densmat)
    if (have_densmat) then
      ASSERT(present(VECUAGR))

      ! Apparently this is only  allocated if N_pc is non-zero. Modern
      ! Fortran  treat  optional arguments  as  not  present when  the
      ! caller passes an unallocated array:
      if (N_pc /= 0 ) then
         ASSERT(present(VECPCGR))
      endif
    else
      ASSERT(present(MATUAGR))

      ! See the comment on present(VECPCGR):
      if (N_pc /= 0 ) then
         ASSERT(present(MATPCGR))
      endif
    endif

    call setif(ISLGR)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    ! also sets LA and LB:
    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! NOT USED: ! allocate global storage of THIS module
    ! NOT USED: call shgi_glob_alloc(NAB,LA,LB)

    allocate(GRSOLV(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    allocate(GRSUM(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    allocate(GAB(NAB,2*LA+1,2*LB+1,9),stat=memstat)
    ASSERT(memstat==0)

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(1,0,0)

    ! S5 overlap integrals:
    call shgi_set_ovrl(LA,LB,1,0,0)  ! (3) in shgi_gr_solv_drv

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    !================================================================
    if (have_densmat) then ! contract integrals with density matrix right away ...
      allocate(PAB(NAB,2*LA+1,2*LB+1),stat=memstat)
      ASSERT(memstat==0)

      ! pack density matrix using CUTOFF that was set in shgi_set_ab():
      do mb=1,2*LB+1
      do ma=1,2*LA+1
        PAB(:,ma,mb) = pack(densmat(:,:,ma,mb),CUTOFF)
      enddo
      enddo

      GRSUM = 0.0_rk
      do up=1,N_points
         z=to_calc_grads%Q(up)
         n_equal=to_calc_grads%n_equal(up)

         do ep=1,n_equal
            x=to_calc_grads%xyz(up,ep,:)
            GRSOLV= 0.0_rk
            call shgi_gr_pc(z,x,GRSOLV)

            ! accumulate to later add to grads wrt A and B:
            GRSUM = GRSUM + GRSOLV

            call shgi_gr_wd_to_abc(GRSOLV,GAB,0)

            ! compute trace(PAB,GAB):
            do i=1,9 ! GAX:GCZ
              VAB(i) = sum( GAB(:,:,:,i) * PAB(:,:,:) )
            enddo
            ! FIXME: we seem to store only VAB(GCX:GCZ), maybe i=GCX,GCZ ???

            ism=to_calc_grads%i_symm_sort(up,ep)
            do ua = 1, n_unique_atoms + N_pc
               if (ua <= n_unique_atoms) then
                  ! NOTE: symmetry type is passed as EC=ISM:
                  call shgi_ve_abc_store(0, 0, 0, 0, UA, ISM, VAB, VECUAGR, -1)
               else
                  ! NOTE: symmetry type is passed as EC=ISM:
                  call shgi_ve_abc_store(0, 0, 0, 0, UA, ISM, VAB, VECPCGR, +1)
               end if
            end do
         enddo
      enddo

      ! store accumulated A and B grads, reuse GAB storage:
      call shgi_gr_wd_to_abc(GRSUM,GAB,0)

      ! compute trace(PAB,GAB):
      do i=1,9 ! GAX:GCZ
        VAB(i) = sum( GAB(:,:,:,i) * PAB(:,:,:) )
      enddo

      call shgi_ve_abc_store(IU2,IE2,IU1,IE1,0,0,VAB,VECUAGR,-1) !ab

      deallocate(PAB,stat=memstat)
      ASSERT(memstat==0)
    !================================================================
    else ! not have densmat, output integral matrices
      GRSUM = 0.0_rk
      do up=1,N_points
         z=to_calc_grads%Q(up)
         n_equal=to_calc_grads%n_equal(up)

         do ep=1,n_equal
            x=to_calc_grads%xyz(up,ep,:)
            GRSOLV= 0.0_rk
            call shgi_gr_pc(z,x,GRSOLV)

            ! accumulate to later add to grads wrt A and B:
            GRSUM = GRSUM + GRSOLV

            call shgi_gr_wd_to_abc(GRSOLV,GAB,0)

            ism=to_calc_grads%i_symm_sort(up,ep)
            do ua = 1, n_unique_atoms + N_pc
               if (ua <= n_unique_atoms) then
                  ! NOTE: symmetry type is passed as EC=ISM:
                  call shgi_gr_abc_store(0, 0, 0, 0, UA, ISM, GAB, MATUAGR, -1)
               else
                  ! NOTE: symmetry type is passed as EC=ISM:
                  call shgi_gr_abc_store(0, 0, 0, 0, UA, ISM, GAB, MATPCGR, +1)
               end if
            end do
         enddo
      enddo

      ! store accumulated A and B grads:
      call shgi_gr_wd_store(IU2, IE2, IU1, IE1, 0, 0, GRSUM, MATUAGR, -1) !ab
    endif
    !================================================================

    deallocate(GRSOLV,stat=memstat)
    ASSERT(memstat==0)

    deallocate(GRSUM,stat=memstat)
    ASSERT(memstat==0)

    deallocate(GAB,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()
    ! NOT USED: ! deallocate global storage of THIS module:
    ! NOT USED: call shgi_glob_free()

    FPP_TIMER_STOP(slgr)
    FPP_TIMER_STOP(pcs)
    FPP_TIMER_STOP(totG)
    FPP_TIMER_STOP(tot)
    call shgi_timing()
  end subroutine shgi_gr_solv_drv

  !****************************************************************
  subroutine shgi_gr_solv_drv_vtn(IU1,IE1,IL1,IU2,IE2,IL2,uas,densmat,VECUAGR,VECMPGR,VECMPTQ)
    use unique_atom_module, only: uat=>unique_atom_type,n_unique_atoms
    use solv_cavity_module, only: to_calc_grads,center2sphere
#ifdef WITH_EFP
    use pointcharge_module, only: pointcharge_N,pca=>pointcharge_array,rcm
    use solv_electrostat_module, only: sphere2center
#endif
    use options_module, only: options_integral_expmax
    use shgi_common, only: CUTOFF
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: shgi_gr_wd_store, shgi_timing, &
                           shgi_gr_abc_store, shgi_gr_wd_to_abc, &
                           shgi_ve_abc_store
    use shgi_pcm, only: shgi_gr_pc
#ifdef WITH_EFP
    use efp_module, only: efp, n_efp
#endif
    implicit none
    type(uat)  , intent(in)      :: uas(:)            ! array of unique atoms, normally all of them
    integer(IK), intent(in)      :: IU1,IE1,IL1,IU2,IE2,IL2
    real(RK)     , intent(in)    :: densmat(:,:,:,:)  ! (NEA,NEB,2*LA+1,2*LB+1), section of the density matrix
    real(RK)     , intent(inout) :: VECUAGR(:)        ! (num of tot-sym grads)            -- grads wrt UAs
    real(RK)     , intent(inout) :: VECMPGR(:)        ! (num of tot-sym grads)            -- grads wrt EFP's
    real(RK)     , intent(inout) :: VECMPTQ(:)        ! (num of tot-sym grads)            -- torqs wrt EFP's
    optional :: VECMPGR,VECMPTQ                       ! -- in "property" runs (forces) we dont need matrices
    ! *** end of interface ***

    integer(IK) :: up,ep,N_points,n_equal,ua,ea,N_pc
    real(RK)    :: x(3),z
    integer(IK) :: memstat
    integer(IK) :: ma,mb,i
    integer(IK) :: c_sphere,t_sphere
#ifdef WITH_EFP
    integer(IK) :: ua1,grp_c1,center_t(2),grp_t,grp_efp
    real(RK)    :: r_rcm(3)
#endif
    real(RK), allocatable :: GRSOLV(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6), W,D
    real(RK), allocatable :: GRSUM(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6), W,D
    real(RK), allocatable :: GAB(:,:,:,:)    ! (NAB,2*LA+1,2*LB+1,9) -- wrt A, B, and C
    real(RK), allocatable :: PAB(:,:,:)      ! (NAB,2*LA+1,2*LB+1)   -- packed densmat
    real(RK)              :: VAB(9)          ! (9) -- wrt A, B, and C, trace(PAB,GAB)
#ifdef WITH_EFP
    real(RK)              :: VAB_t(9)
#endif
    logical               :: have_densmat

    DPRINT  'SHGI: shgi_gr_solv_drv_vtn: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

    FPP_TIMER_START(tot)
    FPP_TIMER_START(totG)
    FPP_TIMER_START(pcs)
    FPP_TIMER_START(slgr)

    have_densmat = .true.
#ifdef WITH_EFP
    if( have_densmat )then
      ASSERT(present(VECMPGR))
      ASSERT(present(VECMPTQ))
    endif
#endif

    call setif(ISLGR)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    ! also sets LA and LB:
    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! NOT USED: ! allocate global storage of THIS module
    ! NOT USED: call shgi_glob_alloc(NAB,LA,LB)

    allocate(GRSOLV(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    allocate(GRSUM(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    allocate(GAB(NAB,2*LA+1,2*LB+1,9),stat=memstat)
    ASSERT(memstat==0)

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(1,0,0)

    ! S5 overlap integrals:
    call shgi_set_ovrl(LA,LB,1,0,0)  ! (3) in shgi_gr_solv_drv

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    N_pc=0
#ifdef WITH_EFP
    if(efp .and. n_efp > 0) N_pc=pointcharge_N
#endif

    N_points=to_calc_grads%n_points

    !================================================================
    if( have_densmat )then ! contract integrals with density matrix right away ...
      allocate(PAB(NAB,2*LA+1,2*LB+1),stat=memstat)
      ASSERT(memstat==0)

      ! pack density matrix using CUTOFF that was set in shgi_set_ab():
      do mb=1,2*LB+1
      do ma=1,2*LA+1
        PAB(:,ma,mb) = pack(densmat(:,:,ma,mb),CUTOFF)
      enddo
      enddo

      GRSUM = 0.0_rk
      do up=1,N_points
#ifdef WITH_EFP
         z=(to_calc_grads%Q(up)+to_calc_grads%Q1(up))*0.5_rk
#else
         z=to_calc_grads%Q(up)
#endif
         n_equal=to_calc_grads%n_equal(up)

         do ep=1,n_equal
            x=to_calc_grads%xyz(up,ep,:)
            t_sphere=to_calc_grads%sphere(up,ep)
#ifdef WITH_EFP
            center_t=sphere2center(t_sphere)
            if(center_t(1) <= n_unique_atoms) then
               grp_t=-1
            else
               grp_t=pca(center_t(1)-N_unique_atoms)%group(center_t(2))
            end if
#endif

            GRSOLV= 0.0_rk
            call shgi_gr_pc(z,x,GRSOLV)

            ! accumulate to later add to grads wrt A and B:
            GRSUM = GRSUM + GRSOLV

            call shgi_gr_wd_to_abc(GRSOLV,GAB,0)

            ! compute trace(PAB,GAB):
            do i=1,9 ! GAX:GCZ
              VAB(i) = sum( GAB(:,:,:,i) * PAB(:,:,:) )
            enddo
            ! FIXME: we seem to store only VAB(GCX:GCZ), maybe i=GCX,GCZ ???

#ifdef WITH_EFP
            grp_efp=0
            do_efp_grad=.false.
#endif
            do ua=1,n_unique_atoms+N_pc
               if(ua <= n_unique_atoms) then
                  do ea=1,uas(ua)%n_equal_atoms
                     c_sphere=center2sphere(ua,ea)

                     if(t_sphere == c_sphere) &
                          call shgi_ve_abc_store(0,0,0,0,UA,EA,VAB,VECUAGR,-1)
                  end do
#ifdef WITH_EFP
               else
                  do_efp_grad=.true.
                  ua1=ua-n_unique_atoms
                  do ea=1,pca(ua1)%n_equal_charges
                     grp_c1=pca(ua1)%group(ea)
                     if(grp_c1 /= grp_efp) then
                        grp_efp=grp_c1
                     else
                        cycle
                     end if

                     if(grp_t == grp_c1) then
                        r_rcm=x-rcm(:,grp_c1)

                        call shgi_ve_abc_store(0,0,0,0,UA1,EA,VAB,VECMPGR,-1)

                        VAB_t=VAB; VAB_t(GCX:GCZ)=vector_product(VAB(GCX:GCZ),r_rcm)
                        call shgi_ve_abc_store(0,0,0,0,UA1,EA,VAB_t,VECMPTQ,+1)
                     end if
                  end do
#endif
               end if
            end do
         enddo
      enddo

#ifdef WITH_EFP
      do_efp_grad=.false.
#endif

      ! store accumulated A and B grads, reuse GAB storage:
      call shgi_gr_wd_to_abc(GRSUM,GAB,0)

      ! compute trace(PAB,GAB):
      do i=1,9 ! GAX:GCZ
        VAB(i) = sum( GAB(:,:,:,i) * PAB(:,:,:) )
      enddo

      call shgi_ve_abc_store(IU2,IE2,IU1,IE1,0,0,VAB,VECUAGR,-1) !ab

      deallocate(PAB,stat=memstat)
      ASSERT(memstat==0)
    endif
    !================================================================

    deallocate(GRSOLV,stat=memstat)
    ASSERT(memstat==0)

    deallocate(GRSUM,stat=memstat)
    ASSERT(memstat==0)

    deallocate(GAB,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()
    ! NOT USED: ! deallocate global storage of THIS module:
    ! NOT USED: call shgi_glob_free()

    FPP_TIMER_STOP(slgr)
    FPP_TIMER_STOP(pcs)
    FPP_TIMER_STOP(totG)
    FPP_TIMER_STOP(tot)
    call shgi_timing()
  end subroutine shgi_gr_solv_drv_vtn

  !*************************************************************
  function vector_product(v1, v2)
    implicit none
    real(RK), intent(in) :: v1(:), v2(:) ! (3)
    real(RK) :: vector_product(3)
    !** End of interface *****************************************

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)
  end function vector_product

  !****************************************************************
  subroutine shgi_gr_Q_solv_drv(IU1,IE1,IL1,IU2,IE2,IL2,uas,densmat,SOLVGR)
    use unique_atom_module, only: uat=>unique_atom_type,n_unique_atoms
    use solv_cavity_module, only: to_calc_grads
    use options_module, only: options_integral_expmax
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_common, only: CUTOFF
    use shgi_utils , only: shgi_gr_wd_store, shgi_timing, &
                           shgi_gr_abc_store, shgi_gr_wd_to_abc, &
                           shgi_ve_abc_store
    use shgi_pcm, only: shgi_gr_pc
    implicit none
    type(uat)  , intent(in)  :: uas(:)            ! array of unique atoms, normally all of them
    integer(IK), intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    real(RK)   , intent(in)  :: densmat(:,:,:,:) ! (NEA,NEB,2*LA+1,2*LB+1), section of the density matrix
    real(RK) , intent(inout) :: SOLVGR(:,:)      ! (n-tot-sym-grads,n-surface-points)
    ! *** end of interface ***

    integer(IK) :: up,ep,N_points,n_equal,ua
    real(RK)    :: x(3)
    integer(IK) :: ism ! symmetry type???
    integer(IK) :: memstat
    integer(IK) :: ma,mb,i
    real(RK), allocatable :: GRSOLV(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6), wrt W,D
    real(RK), allocatable :: GAB(:,:,:,:)    ! (NAB,2*LA+1,2*LB+1,9), wrt A,B,C
    real(RK), allocatable :: PAB(:,:,:)      ! (NAB,2*LA+1,2*LB+1) -- packed densmat
    real(RK)              :: VAB(9)          !                     -- trace(PAB,GAB)

    DPRINT  'SHGI: shgi_gr_solv_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

    FPP_TIMER_START(tot)
    FPP_TIMER_START(totG)
    FPP_TIMER_START(pcs)
    FPP_TIMER_START(slgr)

    call setif(ISLGRT)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    ! also sets LA and LB:
    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! NOT USED: ! allocate global storage of THIS module
    ! NOT USED: call shgi_glob_alloc(NAB,LA,LB)

    allocate(GRSOLV(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    allocate(GAB(NAB,2*LA+1,2*LB+1,9),stat=memstat)
    ASSERT(memstat==0)

    allocate(PAB(NAB,2*LA+1,2*LB+1),stat=memstat)
    ASSERT(memstat==0)

    ! pack density matrix using CUTOFF that was set in shgi_set_ab():
    do mb=1,2*LB+1
    do ma=1,2*LA+1
      PAB(:,ma,mb) = pack(densmat(:,:,ma,mb),CUTOFF)
    enddo
    enddo

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(1,0,0)

    ! S5 overlap integrals:
    call shgi_set_ovrl(LA,LB,1,0,0)  ! (3) in shgi_gr_solv_drv

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    N_points=to_calc_grads%n_points

    do up=1,N_points
       n_equal=to_calc_grads%n_equal(up)

       do ep=1,n_equal
          x=to_calc_grads%xyz(up,ep,:)
          GRSOLV= 0.0_rk
          call shgi_gr_pc(1.0_RK,x,GRSOLV)

          call shgi_gr_wd_to_abc(GRSOLV,GAB,0)

          ! compute trace(PAB,GAB):
          do i=1,9 ! gradients from GAX to GCZ
            VAB(i) = sum(PAB * GAB(:,:,:,i))
          enddo

          ! Unfortunately it has to be done here as we need contribution to
          ! solvation gradients for each surface point (second derivatives case)
          call shgi_ve_abc_store(IU2,IE2,IU1,IE1,0,0,VAB,SOLVGR(:,up),-1)

          ism=to_calc_grads%i_symm_sort(up,ep)
          do ua=1,n_unique_atoms
             ! NOTE: symmetry type is passed as EC=ISM:
             call shgi_ve_abc_store(0,0,0,0,UA,ISM,VAB,SOLVGR(:,up),-1)
          end do
       enddo
    enddo

    deallocate(PAB,stat=memstat)
    ASSERT(memstat==0)

    deallocate(GRSOLV,stat=memstat)
    ASSERT(memstat==0)

    deallocate(GAB,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()
    ! NOT USED: ! deallocate global storage of THIS module:
    ! NOT USED: call shgi_glob_free()

    FPP_TIMER_STOP(slgr)
    FPP_TIMER_STOP(pcs)
    FPP_TIMER_STOP(totG)
    FPP_TIMER_STOP(tot)
    call shgi_timing()
  end subroutine shgi_gr_Q_solv_drv

  !****************************************************************
  !************* DRIVER FOR SECOND DERIVATIVES ********************
  !****************************************************************

  subroutine shgi_sd_solv_drv(IU1,IE1,IL1,IU2,IE2,IL2,uas,densmat,SOLVSD)

    use unique_atom_module, only: uat=>unique_atom_type,n_unique_atoms
    use solv_cavity_module, only: to_calc_grads
    use datatype, only: arrmat4 ! PG storage of grs
    use options_module, only: options_integral_expmax
    use shgi_common, only: CUTOFF
    use shgi_utils , only: shgi_timing, shgi_vd_abc_store, &
         shgi_sd_wd_to_abc, shgi_slvd_abc_store,shgi_slvd_abc_store1, &
         shgi_slvd_abc_store2
    use shgi_pcm, only: shgi_sd_pc
    use shgi_pcm, only: shgi_gr_pc
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)       :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)  , intent(in)       :: uas(:)      ! array of unique atoms, normally all of them
    real(RK)   , intent(in)       :: densmat(:,:,:,:) ! (NEA,NEB,2*LA+1,2*LB+1), section of the density matrix
    real(RK)   , intent(inout)    :: SOLVSD(:,:) ! (num of tot-sym grads^2)
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------
    ! for solvation second derivatives
    real(RK), allocatable :: SDSOLV(:,:,:,:,:)          ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK), allocatable :: GRSOLV(:,:,:,:)            ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: SDSUM(:,:,:,:,:)           ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK), allocatable :: DAB(:,:,:,:,:)             ! (NAB,2*LA+1,2*LB+1,9,9)
    real(RK), allocatable :: PAB(:,:,:)      ! (NAB,2*LA+1,2*LB+1) -- packed densmat
    real(RK)              :: VAB(9,9)        !                     -- trace(PAB,DAB)
    real(RK)              :: GC(9,9)         !                     -- trace(PAB,GAB)
    integer(IK) :: up,ep,N_points,n_equal,ua,ub
    integer(IK) :: memstat
    real(RK)    :: x(3),z
    integer(IK) :: ism ! symmetry type???
    integer(IK) :: ma,mb,i,j

    FPP_TIMER_START(tot)
    FPP_TIMER_START(totD)
    FPP_TIMER_START(pcs)
    FPP_TIMER_START(slsd)
DPRINT  'SHGI: shgi_sd_solv_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! NOT USED: ! allocate global storage of THIS module
    ! NOT USED: call shgi_glob_alloc(NAB,LA,LB)

    allocate(DAB(NAB,2*LA+1,2*LB+1,9,9),stat=memstat)
    ASSERT(memstat==0)

    allocate(SDSUM(NAB,2*LA+1,2*LB+1,6,6),stat=memstat)
    ASSERT(memstat==0)
    SDSUM= 0.0_rk

    allocate(SDSOLV(NAB,2*LA+1,2*LB+1,6,6),stat=memstat)
    ASSERT(memstat==0)

    allocate(GRSOLV(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    allocate(PAB(NAB,2*LA+1,2*LB+1),stat=memstat)
    ASSERT(memstat==0)

    ! pack density matrix using CUTOFF that was set in shgi_set_ab():
    do mb=1,2*LB+1
    do ma=1,2*LA+1
      PAB(:,ma,mb) = pack(densmat(:,:,ma,mb),CUTOFF)
    enddo
    enddo

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(1,1,0)
    ! S5 overlap integrals:
    call shgi_set_ovrl(LA,LB,1,1,0) ! (5) in shgi_sd_solv_drv

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    ! if ever somebody changes the names:
    ASSERT(GCZ-GCX==2)

    N_points=to_calc_grads%n_points
    do up=1,N_points
       z=to_calc_grads%Q(up)
       n_equal=to_calc_grads%n_equal(up)

       do ep=1,n_equal
          x=to_calc_grads%xyz(up,ep,:)

          ! compute second derivatives wrt W and D:
          SDSOLV= 0.0_rk
          call shgi_sd_pc(z,x,SDSOLV)

          ! accumulate second derivatives over all PCs:
          SDSUM=SDSUM+SDSOLV

          ! convert secders of *current* PC to secders wrt A, B, and C:
          call shgi_sd_wd_to_abc(SDSOLV,DAB,0)

          ! contract DAB with density matrix forming trace(PAB,DAB):
          do j=1,9 ! GAX to GCZ
          do i=1,9 ! GAX to GCZ
            VAB(i,j) = sum( PAB(:,:,:) * DAB(:,:,:,i,j) )
          enddo
          enddo

          ! compute *first* derivatives wrt W and D:
          GRSOLV= 0.0_rk
          call shgi_gr_pc(z,x,GRSOLV)

          ! contract first with density matrix forming trace(PAB,GRSOLV),
          ! only need d/dC == - d/dW:
          GC(GCX,GCX) = - sum( PAB(:,:,:) * GRSOLV(:,:,:,GWX) ) ! d/dC = - d/dW
          GC(GCX,GCY) = - sum( PAB(:,:,:) * GRSOLV(:,:,:,GWY) ) ! d/dC = - d/dW
          GC(GCX,GCZ) = - sum( PAB(:,:,:) * GRSOLV(:,:,:,GWZ) ) ! d/dC = - d/dW
          ! FIXME: using 9x9 matrix to store a 3-gradient is confusing

          ism=to_calc_grads%i_symm_sort(up,ep)
          do ua=1,n_unique_atoms
             call shgi_slvd_abc_store(IU2,IE2,IU1,IE1,UA,ISM,VAB,SOLVSD,-1) !AxC, BxC, CxA and CxB
             do ub=1,n_unique_atoms
                call shgi_slvd_abc_store1(ua,ub,ISM,VAB,SOLVSD,-1) !CxC
                call shgi_slvd_abc_store2(ua,ub,ISM,GC ,SOLVSD,-1) !CxC
             end do
          end do
       end do
    end do

    ! now finally store tha accumulated derivatives wrt A and B,
    ! reuse here DAB storage:
    call shgi_sd_wd_to_abc(SDSUM,DAB,0)

    ! contract DAB with density matrix forming trace(PAB,DAB):
    do j=1,9 ! GAX to GCZ
    do i=1,9 ! GAX to GCZ
      VAB(i,j) = sum( PAB(:,:,:) * DAB(:,:,:,i,j) )
    enddo
    enddo

    call shgi_vd_abc_store(IU2,IE2,IU1,IE1,0,0,VAB,SOLVSD,-1)

    deallocate(PAB,stat=memstat)
    ASSERT(memstat==0)

    deallocate(SDSUM,stat=memstat)
    ASSERT(memstat==0)

    deallocate(GRSOLV,stat=memstat)
    ASSERT(memstat==0)

    deallocate(SDSOLV,stat=memstat)
    ASSERT(memstat==0)

    deallocate(DAB,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()
    ! NOT USED: ! deallocate global storage of THIS module:
    ! NOT USED: call shgi_glob_free()

    FPP_TIMER_STOP(slsd)
    FPP_TIMER_STOP(pcs)
    FPP_TIMER_STOP(totD)
    FPP_TIMER_STOP(tot)
    call shgi_timing()
  end subroutine shgi_sd_solv_drv

  !--------------- End of module ----------------------------------
end module shgi_slv

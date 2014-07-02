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
module cpks_xc_resp
  !-------------------------------------------------------------------
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

! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:   &
                   i4_kind &
                 , r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: xc_resp_lda

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  ! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
  ! WARNING!  this is a stripped down test version of the cpksdervs_xc sub  WARNING!
  ! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!

  subroutine xc_resp_lda( job, gb, vl    &
                        , ispin          &
                        , dvdrho         &
                        , dfdrho         &
                        , wts            & 
                        , wts_gr_nuc     &
                        , phi            &
                        , phi_gr_nuc     &
                        , rho_gr_nuc     &
                        , rho_gr_imp     &
                        , eny_sd_imp     &
                        )
    !
    ! Purpose : Used in three regimes: 
    !           a) CPKS RHS,
    !           b) xc-response in CPKS INTERATIONS,
    !           c) FINAL derivatives using CPKS SOLUTION
    !
    use cpksdervs_matrices
    use calc3c_switches
    use error_module, only:MyID
    use cpks_grid_utils
    use datatype, only: arrmat5
    use orbitalstore_module, only: orbital_spin_type
    USE DEBUG, only: NaN
    implicit none
    integer(i4_kind)       , intent(in)    :: job
    integer(i4_kind)       , intent(in)    :: gb
    integer(i4_kind)       , intent(in)    :: vl
    integer(i4_kind)       , intent(in)    :: ispin
    real(r8_kind)          , intent(in)    :: dvdrho(vl,2*ispin-1)
    real(r8_kind)          , intent(in)    :: dfdrho(vl,ispin)
    real(r8_kind)          , intent(in)    :: wts(:)
    type(orbital_spin_type), intent(in)    :: phi(:) ! (n_irr)%o(:vl,eig_dim,partners(irr),ispin)
    type(arrmat5)          , intent(in)    :: phi_gr_nuc(:) ! (n_irr)%m(:vl,eig_dim,partners(irr),n_modes,ispin)
    real(r8_kind)          , intent(in)    :: rho_gr_nuc(:,:,:) ! (:vl,n_modes,ispin) gradient of rho wrt symm nuc
    real(r8_kind)          , intent(out)   :: rho_gr_imp(:,:,:) ! (:vl,n_modes,ispin) impl-grad of rho wrt symm nuc
    real(r8_kind)          , intent(out)   :: eny_sd_imp(:,:,:) ! (n_modes,n_modes,ispin) impl-dervs of energy wrt symm nuc
    ! FIXME: get rid of spin in energy derivatives

    real(r8_kind)          , intent(in)    :: wts_gr_nuc(:,:) ! (:vl,n_modes) gradients of wts wrt symm nuc
    optional :: phi_gr_nuc     ! only for RHS and FINAL
    optional :: rho_gr_nuc     ! only for RHS
    optional :: wts_gr_nuc     ! only for RHS
    optional :: rho_gr_imp     ! only for RHS ( -s1(i,j)/2 ) and FINAL ( u(i,a) )
    optional :: eny_sd_imp     ! only for FINAL
    ! *** end of interface ***

    integer(i4_kind) :: n_irrep             ! := size(phi)
    integer(i4_kind) :: partners(size(phi)) ! (n_irrep)

    real(r8_kind) :: rx(vl,size(cpks,1)) ! response of the density   to perturbation
    real(r8_kind) :: vx(vl,size(cpks,1)) ! response of the potential to tot. sym. pert.
    real(r8_kind) :: v1(vl)              ! response of the potential
    real(r8_kind) :: drv(size(cpks,1),size(cpks,1)) ! matrix of sec ders

    real(r8_kind), allocatable :: H1(:,:)  ! (ni,ni) Ham. gradient occ x occ
    real(r8_kind), allocatable :: Q1(:,:)  ! (ni,na) Ham. gradient occ x vir
    integer(i4_kind) :: istat

    integer(i4_kind) :: i_grad
    integer(i4_kind) :: i_gra1
    integer(i4_kind) :: i_gra2
    integer(i4_kind) :: n_modes ! total number of (symmetric) modes

    integer(i4_kind) :: occ_dim, vir_dim
    integer(i4_kind) :: n_ai ! = occ_dim * vir_dim

    integer(i4_kind) :: i_irr,l,spin

    logical :: q_calc

    ! XXX

    q_calc=.not.cpks_Qxc_calculated

    if(gb==1)then
    print*,MyID,'xc_resp_lda(',job,'): cpks_Qxc_calculated=',cpks_Qxc_calculated
    print*,MyID,'xc_resp_lda(',job,'): q_calc             =',q_calc
    endif

    select case ( job )
    case ( 20 ) ! RHS
      ASSERT(present(phi_gr_nuc))
      ASSERT(present(rho_gr_nuc))
      ASSERT(present(wts_gr_nuc))
      ASSERT(present(rho_gr_imp))
    case ( 24 ) ! FINAL
      ASSERT(present(phi_gr_nuc))
      ASSERT(present(rho_gr_imp))
      ASSERT(present(eny_sd_imp))
    end select

   ASSERT(ispin==1)

   n_modes = size(cpks,1)

   ! for convenient loop ranges:
   n_irrep = size(phi)
   do i_irr=1,n_irrep
     partners(i_irr) = size(phi(i_irr)%o,3)
   enddo

!   n_equal_max=maxval(unique_atoms(:)%n_equal_atoms) ! maximum number of equal_atoms
!   full_vl=machineparameters_veclen ! full vl

    ! note vl is allways smaller or equal then full_vl=machineparameters_veclen
    ! orbs_ob and rho have dimension full_vl

    if_q_calc: if(q_calc) then
    ![[=== Computing RHS Q(a,i) and H1(i,j) for CPKS =======

    if(gb==1)then
    print*,MyID,'xc_resp_lda(',job,'): JOBTYPE = CPKS RHS (req phi_gr_nuc)'
    endif

    spn1: do spin=1,ispin

    FPP_TIMER_START(t_s1_imp_grarho)
    !((=== compute the density response due to overlap changes ===
    rx(:,:) = 0.0_r8_kind
    do i_irr=1,n_irrep
      occ_dim = size(cpks(1,i_irr,spin)%s1,1)
      if( occ_dim == 0 ) cycle ! irrep

      do i_grad=1,size(cpks,1)   
        do l=1,partners(i_irr)
          ! calculate density response 2? * SUM[i,j] { 2 * phi(i) * phi(j) * ( - S1(i,j) / 2 )
          rx(:vl,i_grad) = rx(:vl,i_grad)                         &
                         + rho1( vl, 0                            &
                               , phi(i_irr)%o(:,:,l,spin)         &
                               , - cpks(i_grad,i_irr,spin)%s1 / 2 &
                               )                                  &
                         * (2/ispin)
          ! FIXME: take factros out of args!
        enddo
      enddo
    enddo
    !))===========================================================
    FPP_TIMER_STOP(t_s1_imp_grarho)

    ! output density response due to occ x occ orbital shifts - S1(i,j) / 2:
    rho_gr_imp(:vl,:,spin) = rx(:vl,:)

    ! *** calculate Qai and H1

    irr1: do i_irr=1,n_irrep
!      eig_dim = size(eigvec(i_irr)%m,1)
       occ_dim = size(cpks(1,i_irr,spin)%Qai,1)
       vir_dim = size(cpks(1,i_irr,spin)%Qai,2)
       n_ai    = occ_dim * vir_dim
       if(occ_dim == 0) cycle ! irrep

           ! ** calculate  Q structure  s1 contrib 
           ! ** calculate h1 contibs


  FPP_TIMER_START(t_xc_qh_expl)

  ![[============================================================
  !
  ! Prepare d/dN and integrate over grid
  !
  !       phi(i) * d(w*f)/dN * phi(a)
  ! and
  !       phi(i) * d(w*f)/dN * phi(j)
  !

  ! 1. dRho/dS1: potential response to nuclear displacement or orbital shifts:
  do i_grad=1,size(cpks,1)   
    vx(:,i_grad) = dvdrho(:vl,spin) * rx(:vl,i_grad) &
                 * wts(:vl) / partners(i_irr)
  enddo

  FPP_TIMER_START(t_h1q_dvdrho)

  do i_grad=1,size(cpks,1)
        ! 2. dRho/dN: potential response to nuclear displacement:
        vx(:vl,i_grad) = vx(:vl,i_grad)                                       &
                       + dvdrho(:vl,spin) * ( - rho_gr_nuc(:vl,i_grad,spin) ) &
                       * wts(:vl) / partners(i_irr)
        ! NOTE that rho_gr_nuc = - d/dN!

        ! 3. dW/dN: potential response to integration weight changes:
        vx(:vl,i_grad) = vx(:vl,i_grad)                               &
                       + dfdrho(:vl,spin) * wts_gr_nuc(:vl,i_grad)    &
                                  / partners(i_irr)
        ! FIXME: wts_gr_nuc come already premultiplied with wts!
  enddo

  ! Integrate over grid: 1. Qai( dRho/dS1 ),
  !                      2. Qai( dRho/dN  ),
  !                      3. Qai(   dW/dN  )
  ! (all d/dN -- explicit).
  do i_grad=1,size(cpks,1)
    do l=1,partners(i_irr)
      call gridInt(vl, occ_dim, vx(:,i_grad)    &
                  , phi(i_irr)%o(:,:,l,spin)    &
                  , cpks(i_grad,i_irr,spin)%Qai &
                  )
    enddo
  enddo

  ! Integrate over grid: 1. H1( dRho/dS1 ),
  !                      2. H1( dRho/dN  ),
  !                      3. H1(   dW/dN  )
  ! (all d/dN -- explicit).
  do i_grad=1,size(cpks,1)
    do l=1,partners(i_irr)
      call gridInt(vl, 0, vx(:,i_grad)          &
                  , phi(i_irr)%o(:,:,l,spin)    &
                  , cpks(i_grad,i_irr,spin)%h1  &
                  )
    enddo
  enddo
  !]]============================================================


     FPP_TIMER_STOP(t_h1q_dvdrho)

     ![[=== contributions Phi(i) * (w * df/dr) * d/dN phi(a) ==
     v1(:) = - dfdrho(:vl,spin) * wts(:vl) / partners(i_irr)
     ! FIXME: note minus! Because phi_gr_nuc = - d/dN

     allocate(  H1(occ_dim,occ_dim) &
             ,  Q1(occ_dim,vir_dim) &
             , stat=istat )
     ASSERT(istat==0)

     ! Integrate over grid: phi(i) * (w*f) * d/dN phi(j),
     !                 and  d/dN phi(i) * (w*f) * phi(j):
     ! (d/dN -- explicit).

     do i_grad=1,n_modes

         ! occ x occ:
         H1 = 0.0_r8_kind

         do l=1,partners(i_irr)
           ! do first h1(i,j) := phi(i) (wf) d/dN phi(j)
           call gridInt2(vl, 0, v1(:)                           &
                       , phi(i_irr)%o(:,:,l,spin)               &
                       , phi_gr_nuc(i_irr)%m(:,:,l,i_grad,spin) &
                       , H1                                     &
                       )
         enddo
         ! and then h1(i,j) := h1(i,j) + h1(j,i):
         H1 = H1 + transpose(H1)

         ! occ x vir:
         Q1 = 0.0_r8_kind

         do l=1,partners(i_irr)
           ! do first q1(i,a) := phi(i) (wf) d/dN phi(a)
           call gridInt2(vl, occ_dim, v1(:)                     &
                       , phi(i_irr)%o(:,:,l,spin)               &
                       , phi_gr_nuc(i_irr)%m(:,:,l,i_grad,spin) &
                       , Q1                                     &
                       )

           ! and then q1(i,a) += d/dN phi(i) (wf) phi(a)
           call gridInt2(vl, occ_dim, v1(:)                     &
                       , phi_gr_nuc(i_irr)%m(:,:,l,i_grad,spin) &
                       , phi(i_irr)%o(:,:,l,spin)               &
                       , Q1                                     &
                       )
         enddo

         cpks(i_grad,i_irr,spin)%h1  = cpks(i_grad,i_irr,spin)%h1  + H1
         cpks(i_grad,i_irr,spin)%Qai = cpks(i_grad,i_irr,spin)%Qai + Q1
     enddo ! i_grad

     deallocate( H1 &
               , Q1 &
               , stat=istat )
     ASSERT(istat==0)
     !]]=======================================================

        FPP_TIMER_STOP(t_xc_qh_expl)

    enddo irr1
    enddo spn1
    !]]=== EOF Computing RHS Q(a,i) and H1(i,j) for CPKS ====
  else if_q_calc !----------------------------------
    ![[=== Computing response due to XC kernel in  CPKS iteraltions ===
    !            (and final derivatives in the last iteration)

    ! YYY

    if(gb==1)then
     if( i_cpks.gt.n_cpks )then
       print*,MyID,'xc_resp_lda(',job,'): JOBTYPE = CPKS FINAL (req phi_gr_nuc)'
     else
       print*,MyID,'xc_resp_lda(',job,'): JOBTYPE = CPKS ITERATION'
     endif
    endif

    spn2: do spin=1,ispin

    FPP_TIMER_START(t_xc_bimp_grarho)

    ![[=== compute density response due to orbital shifts U(i,a) ===
    rx(:,:) = 0.0_r8_kind
    do i_irr=1,n_irrep
       occ_dim = size(cpks(1,i_irr,spin)%Qai,1)
       vir_dim = size(cpks(1,i_irr,spin)%Qai,2)
       n_ai    = occ_dim * vir_dim
       if( n_ai == 0 ) cycle

       do l=1,partners(i_irr)
         do i_grad=1,size(cpks,1)

         ! calculates 2? * SUM[i,a] { 2 * phi(i) * phi(a) *  U1(i,a)
         rx(:vl,i_grad) = rx(:vl,i_grad)                    &
                        + rho1( vl, occ_dim                 &
                              , phi(i_irr)%o(:,:,l,spin)    &
                              , cpks(i_grad,i_irr,spin)%HBH &
                              )                             &
                        * (2/ispin)
         enddo
       enddo
    enddo
    !]]=============================================================

    FPP_TIMER_STOP(t_xc_bimp_grarho)

    !=== for final derivatives after CPKS converged ============
    if(i_cpks.gt.n_cpks) then
      ! output density response due to occ x vir orbital shifts U(i,a):
      rho_gr_imp(:vl,:,spin) = rx(:vl,:)

      eny_sd_imp(:,:,spin) = 0.0_r8_kind
    endif

    irr2: do i_irr=1,n_irrep
!      eig_dim = size(eigvec(i_irr)%m,1)
       occ_dim = size(cpks(1,i_irr,spin)%Qai,1)
       vir_dim = size(cpks(1,i_irr,spin)%Qai,2)
       n_ai    = occ_dim * vir_dim
       if ( n_ai == 0 ) cycle irr2 ! next irrep

    ! potential response to nuclear displacement or orbital shifts:
    do i_grad=1,size(rx,2)
      vx(:vl,i_grad) = dvdrho(:vl,spin) * rx(:vl,i_grad) &
                   * wts(:vl) / partners(i_irr)
    enddo

    if( i_cpks <= n_cpks )then
      ! output: kernel response ABi to the current orbital shift
      ! executed many times!

      !=== kernel response integration on the grid ================
      FPP_TIMER_START(t_xc_ab)
      do i_grad=1,size(cpks,1)
        do l=1,partners(i_irr)
           ! Integrate over grid ABi (cpks kernel response):
           call gridInt(vl, occ_dim, vx(:,i_grad)    &
                       , phi(i_irr)%o(:,:,l,spin)    &
                       , cpks(i_grad,i_irr,spin)%ABi &
                       )
        enddo 
      enddo 
      FPP_TIMER_STOP(t_xc_ab)

    else
      ! output: full Fock matrix response: F1 == H1
      ! executed once!

      !=== final response of the fock matrix integration ==========
      do i_grad=1,size(cpks,1)
        do l=1,partners(i_irr)
          ! Integrate over grid h1 (full):
          call gridInt(vl, 0, vx(:,i_grad)          &
                      , phi(i_irr)%o(:,:,l,spin)    &
                      , cpks(i_grad,i_irr,spin)%h1  &
                      )
        enddo 
          ! FIXME: dont compute and invalidate kernel response *after* CPKS:
          cpks(i_grad,i_irr,spin)%ABi = NaN()
      enddo 

      ![[=== for final derivatives after CPKS converged =============
      FPP_TIMER_START(t_dervs_imp)

      !((=== C. compute explicit second derivatives (final) ====

      drv(:,:) = 0.0_r8_kind

      !=== contributions phi(i) * (w * df/dr) * d/dN phi(a) ==
      v1(:) = - dfdrho(:vl,spin) * wts(:vl) !/ partners(i_irr)
      ! FIXME: note minus! Because phi_gr_nuc = - d/dN

      allocate(  H1(occ_dim,occ_dim) &
              ,  Q1(occ_dim,vir_dim) &
              , stat=istat )
      ASSERT(istat==0)

      ! Integrate over grid: phi(i) * (w*f) * d/dN phi(j),
      !                 and  d/dN phi(i) * (w*f) * phi(j):
      ! (d/dN -- explicit).

      do i_gra1=1,n_modes

          ! occ x vir:
          Q1 = 0.0_r8_kind

          do l=1,partners(i_irr)
            ! do first q1(i,a) := phi(i) (wf) d/dN phi(a)
            call gridInt2(vl, occ_dim, v1(:)                     &
                        , phi(i_irr)%o(:,:,l,spin)               &
                        , phi_gr_nuc(i_irr)%m(:,:,l,i_gra1,spin) &
                        , Q1                                     &
                        )

            ! and then q1(i,a) += d/dN phi(i) (wf) phi(a)
            call gridInt2(vl, occ_dim, v1(:)                     &
                        , phi_gr_nuc(i_irr)%m(:,:,l,i_gra1,spin) &
                        , phi(i_irr)%o(:,:,l,spin)               &
                        , Q1                                     &
                        )
          enddo

          ! occ x occ:
          H1 = 0.0_r8_kind

          do l=1,partners(i_irr)
            ! do first h1(i,j) := phi(i) (wf) d/dN phi(j)
            call gridInt2(vl, 0, v1(:)                           &
                        , phi(i_irr)%o(:,:,l,spin)               &
                        , phi_gr_nuc(i_irr)%m(:,:,l,i_gra1,spin) &
                        , H1                                     &
                        )
          enddo
          ! CHOICE: EITHER HERE h1(i,j) := h1(i,j) + h1(j,i):
          H1 = H1 + transpose(H1)

          do i_gra2=1,n_modes
            ! then compute the second derivative of the energy:
            !
            ! drv(X,Y) += 4? * SUM(i,a) UY(i,a)
            !           * ( <      phi(i) | v | d/dX phi(a) >
            !             + < d/dX phi(i) | v |      phi(a) > )
            !
            drv(i_gra1,i_gra2) = drv(i_gra1,i_gra2)                        &
                               - trace2(Q1,cpks(i_gra2,i_irr,spin)%HBH)    &
                               * (4/ispin)
            !
            ! drv(X,Y) += 4? * SUM(i,j) ( - SY(i,j) / 2 )
            !           * ( <      phi(i) | v | d/dX phi(j) >
            !             + < d/dX phi(i) | v |      phi(j) > )
            !
            ! CHOICE: OR A FACTOR 2 HERE:
            drv(i_gra1,i_gra2) = drv(i_gra1,i_gra2)                        &
                               + trace2(H1,cpks(i_gra2,i_irr,spin)%s1) / 2 &
                               * (4/ispin)
            ! NOTE: UX(i,j) == - SX(i,j) / 2 and is symmetric wrt i <-> j
            ! NOTE: remember that phi_gr_nuc = - d/dN phi!
          enddo ! i_gra2
      enddo ! i_gra1

      deallocate( H1 &
                , Q1 &
                , stat=istat )
      ASSERT(istat==0)

      ! output energy derivatives due to the orbital shifts:
      eny_sd_imp(:,:,spin) = eny_sd_imp(:,:,spin) &
                           + drv(:,:)
      !))=======================================================
      FPP_TIMER_STOP(t_dervs_imp)
      !]]=== EOF for final derivatives after CPKS converged =========
    endif

    enddo irr2
    enddo spn2
    !]]=== EOF Computing response due to XC kernel in  CPKS iteraltions ===
 endif if_q_calc

  end subroutine xc_resp_lda

  function trace2(X,Y) result(r)
    !
    ! Compute kind of ``trace of a product''
    !
    !   r = SUM[i,a] X(i,a) * Y(i,a)
    !
    implicit none
    real(r8_kind)   , intent(in) :: X(:,:) ! (ni,na)
    real(r8_kind)   , intent(in) :: Y(:,:) ! (ni,na)
    real(r8_kind)                :: r      ! result
    ! *** end of interface ***

    integer(i4_kind) :: i,a

    ASSERT(size(X,1)==size(Y,1))
    ASSERT(size(X,2)==size(Y,2))

    r = 0.0_r8_kind
    do a=1,size(X,2)
    do i=1,size(X,1)
      r = r + X(i,a) * Y(i,a)
    enddo
    enddo
  end function trace2


  !--------------- End of module -------------------------------------
end module cpks_xc_resp

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
module shgi_sym
  !-------------------------------------------------------------------
  !
  ! Convert  derivatives to  those with  respect to  totally symmetric
  ! modes.
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
  ! Copyright (c) 2006-2008 Alexey Shor
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
# include "def.h"
  use type_module, only: &
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
#ifdef _COMPAC_FORTRAN
  use datatype
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------
!!$  type, public ::  shgi_sym_
!!$  end type shgi_sym_

  !------------ Declaration of constants and variables ---------------
!!$  integer(kind=i4_kind), parameter, public  :: shgi_sym_
!!$  real(kind=r8_kind),    parameter, public  :: shgi_sym_
!!$  logical,               parameter, public  :: shgi_sym_
!!$  character,             parameter, public  :: shgi_sym_
!!$  integer(kind=i4_kind),            public  :: shgi_sym_
!!$  real(kind=r8_kind),               public  :: shgi_sym_
!!$  logical,                          public  :: shgi_sym_
!!$  character,                        public  :: shgi_sym_


  !------------ Interface statements ---------------------------------
!!$  interface shgi_sym_
!!$  end interfaceshgi_sym_
!!$  public shgi_sym_

  !------------ public functions and subroutines ---------------------
  public :: shgi_sym_ve_rot    ! vector version of shgi_sym_gr_rot
  public :: shgi_sym_gr_rot    ! array version of shgi_sym_ve_rot

  public :: shgi_sym_slve_rot  ! vector version of shgi_sym_slgr_rot
  public :: shgi_sym_slgr_rot  ! array version of shgi_sym_slve_rot

  public :: shgi_sym_vd_rot    ! vector version of shgi_sym_sd_rot
  public :: shgi_sym_sd_rot    ! array version of shgi_sym_vd_rot

  public :: shgi_sym_slvd_rot  ! vector version of shgi_sym_slsd_rot

  public :: shgi_sym_fl_rot
  public :: shgi_sym_X_rot

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------
!!$  type
!!$  end type

  !------------ Declaration of constants and variables ---------------
!!$  integer(kind=i4_kind), parameter ::
!!$  real(kind=r8_kind),    parameter ::
!!$  logical,               parameter ::
!!$  character,             parameter ::
!!$  integer(kind=i4_kind),           ::
!!$  real(kind=r8_kind),              ::
!!$  logical,                         ::
!!$  character,                       ::



  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine shgi_sym_X_rot(Xc,UC,EC,GR,SY,NGR)
    !  Purpose: Rotate the X centers grad vector to the
    !           direction of (symmetric) modes (only C center)
    use pointcharge_module, only: unique_pc_index, unique_pc_grad_info
    use point_dqo_module, only: dipoles_grad_index, dipoles_grad_info
    use point_dqo_module, only: qpoles_grad_index, qpoles_grad_info
    use point_dqo_module, only: opoles_grad_index, opoles_grad_info
    use point_dqo_module, only: rep_grad_index, rep_grad_info
    use induced_dipoles_module, only: ind_dip_grad_index,ind_dip_grad_info
    use shgi_cntrl
    implicit none
    integer(IK), intent(in)    :: Xc,UC,EC
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(out)   :: SY(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    integer(IK), intent(out)   :: NGR
    ! *** end of interface ***

    real(RK) :: rot(3,3)

    ASSERT(size(GR,4)==3)
    ASSERT(size(SY,4)==3)

    if(Xc==PC) then
       ngr=unique_pc_index(UC+1)-unique_pc_index(UC)
       rot(:ngr,:)=unique_pc_grad_info(UC)%m(:,:,EC)
    elseif(Xc==PD) then
       ngr=dipoles_grad_index(UC+1)-dipoles_grad_index(UC)
       rot(:ngr,:)=dipoles_grad_info(UC)%m(:,:,EC)
    elseif(Xc==PQ) then
       ngr=qpoles_grad_index(UC+1)-qpoles_grad_index(UC)
       rot(:ngr,:)=qpoles_grad_info(UC)%m(:,:,EC)
    elseif(Xc==PO) then
       ngr=opoles_grad_index(UC+1)-opoles_grad_index(UC)
       rot(:ngr,:)=opoles_grad_info(UC)%m(:,:,EC)
    elseif(Xc==IPD) then
       ngr=ind_dip_grad_index(UC+1)-ind_dip_grad_index(UC)
       rot(:ngr,:)=ind_dip_grad_info(UC)%m(:,:,EC)
    elseif(Xc==RC) then
       ngr=rep_grad_index(UC+1)-rep_grad_index(UC)
       rot(:ngr,:)=rep_grad_info(UC)%m(:,:,EC)
    end if
    SY = mm3(NGR,ROT,GR)
  end subroutine shgi_sym_X_rot

  subroutine shgi_sym_fl_rot(UC,EC,GR,SY,NGR)
    !  Purpose: Rotate the electrostatic field vector to the
    !           direction of (symmetric) modes (only C center)
    use elec_static_field_module, only: surf_points_grad_index
    use elec_static_field_module, only: surf_points_grad_info
    implicit none
    integer(IK), intent(in)    :: UC,EC
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(out)   :: SY(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    integer(IK), intent(out)   :: NGR
    ! *** end of interface ***

    real(RK) :: rot(3,3)

    ASSERT(size(GR,4)==3)
    ASSERT(size(SY,4)==3)

    ngr=surf_points_grad_index(UC+1)-surf_points_grad_index(UC)
    rot(:ngr,:) = surf_points_grad_info(UC)%m(:,:,EC)

    SY = mm3(NGR,ROT,GR)
  end subroutine shgi_sym_fl_rot

  subroutine shgi_sym_slve_rot(UA,ISM,GR,SY,NGR)
    !  Purpose: Rotate the solv gradient to the
    !           direction of (symmetric) modes (only C center)
    use solv_cavity_module, only: to_calc_grads
    use gradient_data_module, only: gradient_index
    use unique_atom_module, only: n_unique_atoms
    use elec_static_field_module, only: surf_points_grad_index
    implicit none
    integer(IK), intent(in)    :: UA,ISM
    real(RK)   , intent(in)    :: GR(:) ! (3)
    real(RK)   , intent(out)   :: SY(:) ! (3)
    integer(IK), intent(out)   :: NGR
    ! *** end of interface ***

    integer(IK) :: j,i

    ASSERT(size(GR)==3)
    ASSERT(size(SY)==3)

    if(ua <= n_unique_atoms) then
       i=ua
       ngr=gradient_index(i+1)-gradient_index(i)
    else
       i=ua-n_unique_atoms
       ngr=surf_points_grad_index(i+1)-surf_points_grad_index(i)
    end if

    do j=1,ngr
       SY(j)=GR(1)*to_calc_grads%dxyz_totsyms(j,UA)%m(1,ISM)+ &
             GR(2)*to_calc_grads%dxyz_totsyms(j,UA)%m(2,ISM)+ &
             GR(3)*to_calc_grads%dxyz_totsyms(j,UA)%m(3,ISM)
    end do
  end subroutine shgi_sym_slve_rot

  subroutine shgi_sym_slgr_rot(UA,ISM,GR,SY,NGR)
    !  Purpose: Rotate the solv gradient to the
    !           direction of (symmetric) modes (only C center)
    use solv_cavity_module, only: to_calc_grads
    use gradient_data_module, only: gradient_index
    use unique_atom_module, only: n_unique_atoms
    use elec_static_field_module, only: surf_points_grad_index
    implicit none
    integer(IK), intent(in)    :: UA,ISM
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(out)   :: SY(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    integer(IK), intent(out)   :: NGR
    ! *** end of interface ***

    integer(IK) :: j,i

    ASSERT(size(GR,4)==3)
    ASSERT(size(SY,4)==3)

    if(ua <= n_unique_atoms) then
       i=ua
       ngr=gradient_index(i+1)-gradient_index(i)
    else
       i=ua-n_unique_atoms
       ngr=surf_points_grad_index(i+1)-surf_points_grad_index(i)
    end if

    do j=1,ngr
       SY(:,:,:,j)=GR(:,:,:,1)*to_calc_grads%dxyz_totsyms(j,UA)%m(1,ISM)+ &
                   GR(:,:,:,2)*to_calc_grads%dxyz_totsyms(j,UA)%m(2,ISM)+ &
                   GR(:,:,:,3)*to_calc_grads%dxyz_totsyms(j,UA)%m(3,ISM)
    end do
  end subroutine shgi_sym_slgr_rot

  subroutine shgi_sym_slvd_rot(U1,E1,U2,E2,IC,GR,SY,NG1,NG2)
    !
    !  Purpose: Rotate the solv second derivatives to the
    !           direction of (symmetric) modes (only CC,CA,CB,AC,BC)
    !           (vector version)
    !
    use solv_cavity_module, only: to_calc_grads
    use gradient_data_module, only: gradient_index
    implicit none
    integer(IK), intent(in)    :: U1, E1
    integer(IK), intent(in)    :: U2, E2
    integer(IK), intent(in)    :: IC
    real(RK)   , intent(in)    :: GR(:,:) ! (3,3)
    real(RK)   , intent(out)   :: SY(:,:) ! (3,3)
    integer(IK), intent(out)   :: NG1
    integer(IK), intent(out)   :: NG2
    ! *** end of interface ***

    integer(IK) :: i,j,ism
    real(RK)    :: ROT1(3,3), ROT2(3,3), BUF(3,3)


    select case( IC )
    case ( 1 )
       ISM=E1
       NG1=gradient_index(U1+1)-gradient_index(U1)
       call gradinfo(U2,E2,ROT2,NG2)
       do i=1,3
          do j=1,NG1
             BUF(j,i)=GR(1,i)*to_calc_grads%dxyz_totsyms(j,U1)%m(1,ism)+ &
                      GR(2,i)*to_calc_grads%dxyz_totsyms(j,U1)%m(2,ism)+ &
                      GR(3,i)*to_calc_grads%dxyz_totsyms(j,U1)%m(3,ism)
          end do
       end do
       do i=1,NG1
          SY(i,:) = ro3(NG2,ROT2,BUF(i,:))
       end do
    case ( 2 )
       ISM=E2
       call gradinfo(U1,E1,ROT1,NG1)
       NG2=gradient_index(U2+1)-gradient_index(U2)
       do i=1,3
          BUF(:,i) = ro3(NG1,ROT1,GR(:,i))
       end do
       do i=1,NG1
          do j=1,NG2
             SY(i,j)=BUF(i,1)*to_calc_grads%dxyz_totsyms(j,U2)%m(1,ism)+ &
                     BUF(i,2)*to_calc_grads%dxyz_totsyms(j,U2)%m(2,ism)+ &
                     BUF(i,3)*to_calc_grads%dxyz_totsyms(j,U2)%m(3,ism)
          end do
       end do
    case ( 3 )
       ISM=E1 !or E2
       NG1=gradient_index(U1+1)-gradient_index(U1)
       NG2=gradient_index(U2+1)-gradient_index(U2)
       do i=1,3
          do j=1,NG1
             BUF(j,i)=GR(1,i)*to_calc_grads%dxyz_totsyms(j,U1)%m(1,ism)+ &
                      GR(2,i)*to_calc_grads%dxyz_totsyms(j,U1)%m(2,ism)+ &
                      GR(3,i)*to_calc_grads%dxyz_totsyms(j,U1)%m(3,ism)
          end do
       end do
       do i=1,NG1
          do j=1,NG2
             SY(i,j)=BUF(i,1)*to_calc_grads%dxyz_totsyms(j,U2)%m(1,ism)+ &
                     BUF(i,2)*to_calc_grads%dxyz_totsyms(j,U2)%m(2,ism)+ &
                     BUF(i,3)*to_calc_grads%dxyz_totsyms(j,U2)%m(3,ism)
          end do
       end do
    case ( 4 )
       ISM=E1 !or E2
       NG1=gradient_index(U1+1)-gradient_index(U1)
       NG2=gradient_index(U2+1)-gradient_index(U2)

       do i=1,NG1
          do j=1,NG2
             SY(i,j)=GR(1,1)*to_calc_grads%d2xyz_totsyms(i,U1,j,U2)%m(1,ism)+ &
                     GR(1,2)*to_calc_grads%d2xyz_totsyms(i,U1,j,U2)%m(2,ism)+ &
                     GR(1,3)*to_calc_grads%d2xyz_totsyms(i,U1,j,U2)%m(3,ism)
          end do
       end do
    case default
       ABORT('no such case')
    end select
  end subroutine shgi_sym_slvd_rot

  subroutine shgi_sym_ve_rot(UA,EA,GR,SY,NGR)
    !
    !  Purpose: Rotate the gradient vector to the
    !           direction of (symmetric) modes
    !
    implicit none
    integer(IK), intent(in)    :: UA, EA
    real(RK)   , intent(in)    :: GR(:) ! (3)
    real(RK)   , intent(out)   :: SY(:) ! (3)
    integer(IK), intent(out)   :: NGR
    ! *** end of interface ***

    real(RK)    :: ROT(3,3)
    integer(IK) :: j

    ASSERT(size(GR)==3)
    ASSERT(size(SY)==3)

    call gradinfo(UA,EA,ROT,NGR)

    do j=1,NGR
      SY(j) = ROT(j,1)*GR(1) &
            + ROT(j,2)*GR(2) &
            + ROT(j,3)*GR(3)
    enddo
    ! SY(NGR+1:) is undefined in case NGR < 3!
  end subroutine shgi_sym_ve_rot

  subroutine shgi_sym_gr_rot(UA,EA,GR,SY,NGR)
    !
    !  Purpose: Rotate the gradient matrix to the
    !           direction of (symmetric) modes
    !
    !  Note: if possible, first contract the matrix
    !        with the density matrix to form a vector
    !        and then rotate...
    !
    implicit none
    integer(IK), intent(in)    :: UA, EA
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(out)   :: SY(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    integer(IK), intent(out)   :: NGR
    ! *** end of interface ***

    real(RK)    :: ROT(3,3)

    ASSERT(size(GR,4)==3)
    ASSERT(size(SY,4)==3)

    call gradinfo(UA,EA,ROT,NGR)

    SY = mm3(NGR,ROT,GR)
  end subroutine shgi_sym_gr_rot

  function ro3(N,A,B) result(C)
    !
    ! Vector version of mm3
    !
    implicit none
    real(RK)   , intent(in)    :: A(3,3)
    real(RK)   , intent(in)    :: B(:)
    real(RK)                   :: C(3)
    integer(IK), intent(in)    :: N
    ! *** end of interface ***

    integer(IK) :: j

    do j=1,N
      C(j) = A(j,1)*B(1) &
           + A(j,2)*B(2) &
           + A(j,3)*B(3)
    enddo
  end function ro3

  function mm3(N,A,B) result(C)
    !
    ! Array version of ro3
    !
!   use debug, only: NaN
    implicit none
    real(RK)   , intent(in)    :: A(3,3)
    real(RK)   , intent(in)    :: B(:,:,:,:)
    real(RK)                   :: C(size(B,1),size(B,2),size(B,3),3)
    integer(IK), intent(in)    :: N
    ! *** end of interface ***

    integer(IK) :: j

    do j=1,N
      C(:,:,:,j) = A(j,1)*B(:,:,:,1) &
                 + A(j,2)*B(:,:,:,2) &
                 + A(j,3)*B(:,:,:,3)
    enddo
!ifndef _COMPAC_FORTRAN
!   ! FIXME: for debuging mark the values:
!   do j=N+1,3
!     C(:,:,:,j) = NaN()
!   enddo
!endif
  end function mm3

  subroutine shgi_sym_vd_rot(U1,E1,U2,E2,GR,SY,NG1,NG2)
    !
    !  Purpose: Rotate the gradient to the
    !           direction of (symmetric) modes
    !           (vector version)
    !
    implicit none
    integer(IK), intent(in)    :: U1, E1
    integer(IK), intent(in)    :: U2, E2
    real(RK)   , intent(in)    :: GR(:,:) ! (3,3)
    real(RK)   , intent(out)   :: SY(:,:) ! (3,3)
    integer(IK), intent(out)   :: NG1
    integer(IK), intent(out)   :: NG2
    ! *** end of interface ***

    integer(IK) :: i
    real(RK)    :: ROT1(3,3), ROT2(3,3)

    ASSERT(size(GR,1)==3)
    ASSERT(size(SY,1)==3)
    ASSERT(size(GR,2)==3)
    ASSERT(size(SY,2)==3)

    call gradinfo(U1,E1,ROT1,NG1)
    call gradinfo(U2,E2,ROT2,NG2)

    do i=1,3
      SY(:,i) = ro3(NG1,ROT1,GR(:,i))
    enddo
    do i=1,NG1
      SY(i,:) = ro3(NG2,ROT2,SY(i,:))
    enddo
  end subroutine shgi_sym_vd_rot

  subroutine shgi_sym_sd_rot(U1,E1,U2,E2,GR,SY,NG1,NG2)
    !
    !  Purpose: Rotate the gradient to the
    !           direction of (symmetric) modes
    !           (array version)
    !
    implicit none
    integer(IK), intent(in)    :: U1, E1
    integer(IK), intent(in)    :: U2, E2
    real(RK)   , intent(in)    :: GR(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    real(RK)   , intent(out)   :: SY(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    integer(IK), intent(out)   :: NG1
    integer(IK), intent(out)   :: NG2
    ! *** end of interface ***

    integer(IK) :: i
    real(RK)    :: ROT1(3,3), ROT2(3,3)

    ASSERT(size(GR,4)==3)
    ASSERT(size(SY,4)==3)
    ASSERT(size(GR,5)==3)
    ASSERT(size(SY,5)==3)

    call gradinfo(U1,E1,ROT1,NG1)
    call gradinfo(U2,E2,ROT2,NG2)

    do i=1,3
      SY(:,:,:,:,i) = mm3(NG1,ROT1,GR(:,:,:,:,i))
    enddo
    do i=1,NG1
      SY(:,:,:,i,:) = mm3(NG2,ROT2,SY(:,:,:,i,:))
    enddo
  end subroutine shgi_sym_sd_rot

  subroutine gradinfo(UA,EA,rot,NGR)
    !  Purpose: Rotate the gradient to the
    !           direction of (symmetric) modes
    use unique_atom_module, only: unique_atom_grad_info, n_unique_atoms
#ifdef WITH_EFP
    use pointcharge_module, only: unique_pc_grad_info, pointcharge_N
    use shgi_cntrl, only: do_efp_grad
#endif
    implicit none
    integer(IK), intent(in)  :: UA, EA
    real(RK)   , intent(out) :: rot(3,3)
    integer(IK), intent(out) :: NGR
    ! *** end of interface ***

    integer(IK) :: n_ua
    integer(IK) :: sh2(2)
#ifdef WITH_EFP
    real(RK), pointer :: grad_info(:,:)
#endif

#ifdef WITH_EFP
    if(do_efp_grad) then
       n_ua = size(unique_pc_grad_info)
       ASSERT(n_ua==pointcharge_N)
       grad_info=>unique_pc_grad_info(UA)%m(:,:,EA)
    else
       n_ua = size(unique_atom_grad_info)
       ! moving atoms not yet implemented:
       ASSERT(n_ua==n_unique_atoms)
       grad_info=>unique_atom_grad_info(UA)%m(:,:,EA)
    end if

    sh2 = shape(grad_info)
    ASSERT(sh2(2)==3) ! x,y,z
    ASSERT(sh2(1)<=3) ! tot sym modes

    NGR = sh2(1)
    rot = 0.0_rk
    rot(:NGR,:) = grad_info
#else
    n_ua = size(unique_atom_grad_info)
    ! moving atoms not yet implemented:
    ASSERT(n_ua==n_unique_atoms)

    sh2 = shape( unique_atom_grad_info(UA)%m(:,:,EA) )
    ASSERT(sh2(2)==3) ! x,y,z
    ASSERT(sh2(1)<=3) ! tot sym modes

    NGR = sh2(1)
    rot = 0.0_rk
    rot(:NGR,:) = unique_atom_grad_info(UA)%m(:,:,EA)
#endif
  end subroutine gradinfo

  !--------------- End of module -------------------------------------
end module shgi_sym

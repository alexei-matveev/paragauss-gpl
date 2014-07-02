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
module gga_response_module
  !-------------------------------------------------------------------
  !
  !  Purpose: Contains routines for evaluating PW91 and PBE
  !   exchange-correlation functionals:
  !
  !  Subroutine *x_calc()  -> GGA exchange    of PW91 & PBE (moved to PW91)
  !  Subroutine *c_calc()  -> GGA correlation of PW91 & PBE
  !
  !  Module called by: response_module
  !
  !  References:
  !   1. Based on "svwn_module.f90" module of HH, UB.
  !   2. B.G. Johnson, P.M.W. Gill and. J.A. Pople,
  !      "The performance of a family of density functional methods",
  !      J.Chem.Phys. Vol.98 No.7, April 1993 pp. 5612.
  !      (Erratum: JCP Vol.1 No.10, Nov. 15, pp. 9202)
  !
  !  Author: SB
  !  Date: 5/04
  !
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification
  ! Author: SB
  ! Date:   17/04
  ! Description: Added PBE
  ! Reference:
  !  1. Based on J.P. Perdew et al. "GGA made simple" Phys. Rev. Lett.
  !                                 77(18) P.3865 1996
  !
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: SB
  ! Date:   28.6.2006
  ! Description: rebuild all and change the interface
  !
  !! VN variant (will be implemented and used everywhere)
  !! NOTE for second derivatives packing:
  !!              1    2    3    4    5    6
  !! dfdg        AA   BB   AB
  !! dfdndn      AA   BB   AB
  !! dfdndg     AAA  BBB  BAA  ABB  AAB  BAB
  !! dfdgdg    AAAA BBBB AABB AAAB BBAB ABAB
  !
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification
  ! Author: TS
  ! Date:   11.09.09
  ! Description: Changed phi to phi + eps to avoid singularities
  !              in derivatives of phi for cases where:
  !                      rho = (rho_alpha,0)
  !
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
  use type_module ! type specification parameters
  use constants
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ public functions and subroutines ---------------------
  public gga_correlation

  !-------------------------------------------------------------------
  !! BLOCK FOR PBE FINISHED

  integer(i4_kind), parameter, public :: &
       C_PBE  = 3, &
       C_PW91 = 4


  integer(i4_kind), parameter :: &
       DA   = 1, &
       DB   = 2, &
       DG   = 3, &
       DGDA = 1, &
       DGDB = 2, &
       DADA = 1, &
       DBDB = 2, &
       DADB = 3
  !!       DBDA = 4

  real(kind=r8_kind),parameter :: eps = 1.0E-20_r8_kind


  integer(i4_kind), parameter :: & !! Hln, dHln, d2Hln
       DT      = 1,&
       DX      = 2,&
       DTDT    = 1,&
       DXDX    = 2,&
       DTDX    = 3

  integer(i4_kind), parameter :: & !! A, dA, d2A
       DP      = 1,&
       DE      = 2,&
       DPDP    = 1,&
       DEDE    = 2,&
       DEDP    = 3

  integer(i4_kind), parameter :: &
       DPHI    = 1,&
       DT2     = 2,&
       DEC     = 3,&
       DRS     = 4 !! for PW91

  !! pw91 parameters
  real(kind=r8_kind),parameter :: &
       alpha = 0.09_r8_kind, &
       nu    = 15.755920349483144658863574030524_r8_kind , &
       !!(four*four/pi)*(three*pi*pi)**(one/three), &
       Cc0   = 0.004235_r8_kind, &
       Cx    = - 0.001667_r8_kind

  !------------ Subroutines ------------------------------------------
contains

  subroutine gga_correlation(C_KIND,ispin,id,rho,gmma,vlen,&
       f,Ec,&
       dfdn,dfdg,dEcdn,&
       dfdndn,dfdndg,dfdgdg,dEcdndn)
    ! Purpose: GGA exchange of PBE
    !
    implicit none

    integer(kind=i4_kind), intent(in)    :: C_KIND
    integer(kind=i4_kind), intent(in)    :: ispin
    integer(kind=i4_kind), intent(in)    :: id   !! index of derivatives
    integer(kind=i4_kind), intent(in)    :: vlen
    real(kind=r8_kind),    intent(in)    :: rho(:,:), gmma(:,:)
    real(kind=r8_kind),    intent(inout) :: f(:),dfdn(:,:),dfdg(:,:),dfdndn(:,:),dfdndg(:,:),dfdgdg(:,:)
    real(kind=r8_kind),    intent(in)    :: Ec(:),dEcdn(:,:),dEcdndn(:,:) !! LDA from outside

    optional                             :: dfdn,dfdg,dEcdn
    optional                             :: dfdndn,dfdndg,dfdgdg,dEcdndn

    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind),dimension(vlen,2)     :: nn
    real(kind=r8_kind),dimension(vlen)       :: gm

    integer(i4_kind), parameter ::&
         & AA = 1, &
         & BB = 2, &
         & AB = 3
    !!         & BA = 4

    integer(i4_kind), parameter ::&
         & AAA = 1, &
         & BBB = 2, &
         & BAA = 3, &
         & ABB = 4, &
         & AAB = 5, &
         & BAB = 6

    integer(i4_kind), parameter ::&
         & AAAA = 1, &
         & BBBB = 2, &
         & AABB = 3, &
         & AAAB = 4, &
         & BBAB = 5, &
         & ABAB = 6

!!$    integer(i4_kind) :: i

    !! TDDFT: need to be half of rho

    if (ispin==1) then
       nn(:, 1) = rho(:vlen, 1) / 2
       nn(:, 2) = nn(:, 1)
       gm = gmma(:vlen, 1)
    else
       nn(:, 1) = rho(:vlen, 1)
       nn(:, 2) = rho(:vlen, 2)
       gm = gmma(:vlen, 1) + gmma(:vlen, 2) + TWO * gmma(:vlen, 3)
       where(gm<eps) gm = eps
    endif

    !! ispin       result
    !!  1       dfdn(1:vlen,DA)
    !!  2       dfdn(1:vlen,DB)
    !! 1.or.2   dfdg(1:vlen,DG)
    !!  1       dfdndn(1:vlen,DADA), dfdndn(1:vlen,DADB)
    !!  2       dfdndn(1:vlen,DBDB), dfdndn(1:vlen,DBDA)
    !!  1       dfdndg(1:vlen,DADG)
    !!  2       dfdndg(1:vlen,DBDG)
    !! 1.or.2   dfdgdg(1:vlen,DGDG)

    select case(id)
    case (1)
       call pbe_or_pw91(C_KIND,id,nn,gm,vlen,&
            f,Ec, &
            dfdn,dfdg(:,1),dEcdn)
    case (2)
       call pbe_or_pw91(C_KIND,id,nn,gm,vlen,&
            f,Ec, &
            dfdn,dfdg(:,1),dEcdn,&
            dfdndn,dfdndg(1:vlen,1:2),dfdgdg(1:vlen,1),dEcdndn)
    end select
    !! REBUILDING

    if (ispin==1) return

    !      dfdg  (1:vlen,AA)  = dfdg(1:vlen,1)
    dfdg  (1:vlen,BB)  = dfdg(1:vlen,AA)
    dfdg  (1:vlen,AB)  = dfdg(1:vlen,AA) * TWO

    if (id == 1 ) return

    !      dfdndg(1:vlen,AAA) =       dfdndg(1:vlen,AAA)
    dfdndg(1:vlen,ABB) =       dfdndg(1:vlen,AAA)
    dfdndg(1:vlen,AAB) = TWO * dfdndg(1:vlen,AAA)
    !      dfdndg(1:vlen,BBB) =       dfdndg(1:vlen,BBB)
    dfdndg(1:vlen,BAA) =       dfdndg(1:vlen,BBB)
    dfdndg(1:vlen,BAB) = TWO * dfdndg(1:vlen,BBB)

    !      dfdgdg(1:vlen,AAAA) =        dfdgdg(1:vlen,AAAA)
    dfdgdg(1:vlen,BBBB) =        dfdgdg(1:vlen,AAAA)
    dfdgdg(1:vlen,AABB) =        dfdgdg(1:vlen,AAAA)
    dfdgdg(1:vlen,AAAB) =  TWO * dfdgdg(1:vlen,AAAA)
    dfdgdg(1:vlen,BBAB) =  TWO * dfdgdg(1:vlen,AAAA)
    dfdgdg(1:vlen,ABAB) = FOUR * dfdgdg(1:vlen,AAAA)

  end subroutine gga_correlation
  !********************************************************************

  subroutine pbe_or_pw91(C_KIND,id,rho,gm,vlen,&
       fc,Ec,&
       dfcdn,dfcdg,dEcdn,&
       dfc_drhodrho,dfc_drhodgamma,dfc_dgammadgamma,dEcdndn)
    ! Purpose: GGA correlation for PBE
    !
    use constants
    use pw_ldac_module
    implicit none

    integer(kind=i4_kind), intent(in)    :: C_KIND
    integer(kind=i4_kind), intent(in)    :: id !! index of derivatives
    integer(kind=i4_kind), intent(in)    :: vlen !!, ispin
    real(kind=r8_kind),    intent(in)    :: rho(:,:),gm(:)
    real(kind=r8_kind),    intent(inout) :: fc(:),dfcdn(:,:)
    real(kind=r8_kind),    intent(inout) :: dfcdg(:),dfc_drhodrho(:,:),dfc_drhodgamma(:,:),dfc_dgammadgamma(:)
    real(kind=r8_kind),    intent(in)    :: Ec(:),dEcdn(:,:),dEcdndn(:,:)

    optional                             :: dfcdn,dfcdg,dEcdn
    optional                             :: dfc_drhodrho,dfc_drhodgamma,dfc_dgammadgamma,dEcdndn

    !------------ Declaration of local variables ---------------------
    real(kind=r8_kind),dimension(vlen)                   :: H,dHdgdg
    real(kind=r8_kind),dimension(vlen,2)                 :: dHdndg
    real(kind=r8_kind),dimension(vlen,3)                 :: dHdn, dHdndn
    real(kind=r8_kind),dimension(vlen)                   :: nn
    real(kind=r8_kind),dimension(vlen)                   :: z,dzda,dzdada,dzdb,dzdadb,dzdbdb
    real(kind=r8_kind),dimension(vlen)                   :: phi,dtdphi,dtdphidphi,dtdn,dtdphidg
    real(kind=r8_kind),dimension(vlen)                   :: t,dphidz,dphidzdz
    real(kind=r8_kind),dimension(vlen)                   :: dtdndn,dtdndphi,dtdg,dtdndg

    real(kind=r8_kind),dimension(vlen,DPHI:DRS)          :: dH
    real(kind=r8_kind),dimension(vlen,DPHI:DRS,DPHI:DRS) :: d2H
    real(kind=r8_kind),dimension(vlen,DPHI:DRS,3  )      :: dXdx
    real(kind=r8_kind),dimension(vlen,DPHI:DRS,3,3)      :: d2Xdx2

    real(kind=r8_kind),dimension(vlen)                   ::  Flda
    real(kind=r8_kind),dimension(vlen,2)                 :: dFldadn
    real(kind=r8_kind),dimension(vlen,3)                 :: dFldadndn

    real(kind=r8_kind),dimension(vlen)                   :: rs, drsdn, drsdndn
    integer(kind=i4_kind)                                :: END

    integer(kind=i4_kind)                                :: i,j

    real(kind=r8_kind) :: t2_fact

    real(kind=r8_kind), parameter :: &
         c1d3  = one/three, &
         c1d9  = one/nine,  &
         c2d3  = two/three, &
         c4d3  = four/three,&
         c5d3  = five/three,&
         c7d6  = seven/six

    nn = rho(1:vlen,1)+rho(1:vlen,2)

    where(abs(nn)<eps) nn = eps

    if (C_KIND == C_PW91) then
       rs      = 1.91916_r8_kind / (three * pi * pi * nn) ** c1d3
       drsdn   = - c1d3 * rs     / nn
       drsdndn = - c4d3 * drsdn  / nn
    end if

    !! dzeta derivatives

    z = (rho(1:vlen,1)-rho(1:vlen,2))/nn

    if (id > 0) then
       dzda   =    (one-z)/nn
       dzdb   =   -(one+z)/nn
    end if

    if (id==2) then
       dzdada =  -two*dzda/nn
       dzdadb =   two*   z/(nn*nn)
       dzdbdb =  -two*dzdb/nn
    end if

    !! phi derivatives
    phi      =  half*((one+z) **   c2d3  + (one-z)**  c2d3)
    if (id > 0) then
       dphidz   =  c1d3*((one+z+eps) ** (-c1d3) - (one-z+eps)**(-c1d3))
    end if

    if (id == 2) then
       dphidzdz = -c1d9*((one+z+eps) ** (-c4d3) + (one-z+eps)**(-c4d3))
    end if

    !! t2 derivatives
    t2_fact    = pi / (four*four*(three*pi**2)**c1d3)
    t          = t2_fact*gm/((phi*nn**c7d6)**two)
    if (id > 0 ) then
       dtdphi     = - two   *      t/phi
       dtdn       = -(seven/three) *    t/nn
       dtdg       = t2_fact   /((phi*nn**c7d6)**two)
    end if
    if (id==2) then
       dtdphidphi = - three * dtdphi/phi
       dtdndn     = -(  ten/three) * dtdn/nn
       dtdndg     = -(seven/three) * dtdg/nn
       dtdphidg   = -two * dtdg/phi
       dtdndphi   = -two * dtdn/phi
    end if

    Flda(1:vlen) = Ec(1:vlen)/nn(1:vlen)

    if (id > 0 ) then
       dFldadn(1:vlen,DA)     = dEcdn(1:vlen,DA)    /nn - Flda(1:vlen)/nn
       dFldadn(1:vlen,DB)     = dEcdn(1:vlen,DB)    /nn - Flda(1:vlen)/nn
    end if

    if (id==2) then
       dFldadndn(1:vlen,DADA) = dEcdndn(1:vlen,DADA)/nn - two * dFldadn(1:vlen,DA)/nn
       dFldadndn(1:vlen,DADB) = dEcdndn(1:vlen,DADB)/nn - (dFldadn(1:vlen,DA) + dFldadn(1:vlen,DB))/nn
       dFldadndn(1:vlen,DBDB) = dEcdndn(1:vlen,DBDB)/nn - two * dFldadn(1:vlen,DB)/nn
    end if

    H   = zero
    dH  = zero
    d2H = zero

    select case (C_KIND)
    case (C_PBE)
       call H_pbe (id,   phi,t,Flda,vlen, H, dH, d2H)
    case (C_PW91)
       call H_pw91(id,rs,phi,t,Flda,vlen, H, dH, d2H)
    case default
       STOP
    end select

    ! derivatives (first, and second) of the arguments
    ! Ec, phi, t
    ! first : da, db
    ! second: dada, dadb, dbdb
    ! dXdx  (1:3,1:2)
    ! d2Xdx2(1:3,1:4)
    if (id > 0 ) then
       ! Changed from dphidz*dzda and dphidz*dzdb to allow a proper calling from tpss module
       dXdx(1:vlen,DPHI,DA)      = c1d3 * two**c2d3 / nn**c5d3 * (rho(1:vlen,2) / (rho(1:vlen,1)+eps)**c1d3 - rho(1:vlen,2)**c2d3)
       dXdx(1:vlen,DPHI,DB)      = c1d3 * two**c2d3 / nn**c5d3 * (rho(1:vlen,1) / (rho(1:vlen,2)+eps)**c1d3 - rho(1:vlen,1)**c2d3)
       dXdx(1:vlen,DPHI,DG)      = zero

       dXdx(1:vlen,DT2,  DA)      = dtdn + dtdphi * dXdx(1:vlen,DPHI,DA)
       dXdx(1:vlen,DT2,  DB)      = dtdn + dtdphi * dXdx(1:vlen,DPHI,DB)
       dXdx(1:vlen,DT2,  DG)      = dtdg

       dXdx(1:vlen,DEC, DA)      = dFldadn(1:vlen,DA)
       dXdx(1:vlen,DEC, DB)      = dFldadn(1:vlen,DB)
       dXdx(1:vlen,DEC, DG)      = zero

       if (C_KIND==C_PW91) then
          dXdx(1:vlen,DRS, DA)   = drsdn
          dXdx(1:vlen,DRS, DB)   = drsdn
          dXdx(1:vlen,DRS, DG)   =  zero
       end if

    end if

    if (id==2) then

       d2Xdx2(1:vlen,DPHI,DA,DA) = dphidzdz*dzda*dzda + dphidz*dzdada
       d2Xdx2(1:vlen,DPHI,DA,DB) = dphidzdz*dzda*dzdb + dphidz*dzdadb
       d2Xdx2(1:vlen,DPHI,DB,DB) = dphidzdz*dzdb*dzdb + dphidz*dzdbdb
       !!    d2Xdx2(1:vlen,DPHI,DB,DA) = dphidzdz*dzdb*dzda + dphidz*dzdadb
       d2Xdx2(1:vlen,DPHI,DG,DA) = zero
       d2Xdx2(1:vlen,DPHI,DG,DB) = zero
       d2Xdx2(1:vlen,DPHI,DG,DG) = zero

       d2Xdx2(1:vlen,DT2,DA,DA)   = &
            dtdphidphi*dphidz*dzda * dphidz * dzda + &
            TWO*dtdndphi * dphidz * dzda           + &
            dtdphi * dphidzdz * (dzda**TWO)        + &
            dtdphi * dphidz * dzdada               + &
            dtdndn

       d2Xdx2(1:vlen,DT2,DA,DB)   = dtdndn + dtdndphi * dphidz * dzdb &
            + (dtdndphi + dtdphidphi*dphidz*dzdb) * dphidz * dzda &
            + dtdphi * dphidzdz * dzdb * dzda + dtdphi * dphidz * dzdadb

       d2Xdx2(1:vlen,DT2,DB,DB)   = &
            dtdphidphi*dphidz*dzdb * dphidz * dzdb + &
            TWO*dtdndphi * dphidz * dzdb           + &
            dtdphi * dphidzdz * (dzdb**TWO)        + &
            dtdphi * dphidz * dzdbdb               + &
            dtdndn

       !!    d2Xdx2(1:vlen,DT2,DB,DA)   = dtdndn + dtdndphi * dphidz * dzdb &
       !!         + (dtdndphi + dtdphidphi*dphidz*dzdb) * dphidz * dzda &
       !!         + dtdphi * dphidzdz * dzdb * dzda + dtdphi * dphidz * dzdadb

       d2Xdx2(1:vlen,DT2,DG,DA)   = dtdndg+dtdphidg*dphidz*dzda
       d2Xdx2(1:vlen,DT2,DG,DB)   = dtdndg+dtdphidg*dphidz*dzdb
       d2Xdx2(1:vlen,DT2,DG,DG)   = zero

       d2Xdx2(1:vlen,DEC,DA,DA)  = dFldadndn(1:vlen,DADA)
       d2Xdx2(1:vlen,DEC,DA,DB)  = dFldadndn(1:vlen,DADB)
       d2Xdx2(1:vlen,DEC,DB,DB)  = dFldadndn(1:vlen,DBDB)
       !!    d2Xdx2(1:vlen,DEC,DB,DA)  = dFldadndn(1:vlen,DBDA)
       d2Xdx2(1:vlen,DEC,DG,DA)  = zero
       d2Xdx2(1:vlen,DEC,DG,DB)  = zero
       d2Xdx2(1:vlen,DEC,DG,DG)  = zero

       if (C_KIND==C_PW91) then
          d2Xdx2(1:vlen,DRS, DA,DA) = drsdndn
          d2Xdx2(1:vlen,DRS, DA,DB) = drsdndn
          d2Xdx2(1:vlen,DRS, DB,DB) = drsdndn

!!          d2Xdx2(1:vlen,DRS, DA,DG) = zero
!!          d2Xdx2(1:vlen,DRS, DB,DA) = drsdndn
!!          d2Xdx2(1:vlen,DRS, DB,DG) = zero

          d2Xdx2(1:vlen,DRS, DG,DA) = zero
          d2Xdx2(1:vlen,DRS, DG,DB) = zero
          d2Xdx2(1:vlen,DRS, DG,DG) = zero
       end if

    end if

    dHdn   = zero
    dHdndn = zero
    dHdndg = zero
    dHdgdg = zero

    fc(1:vlen) = fc(1:vlen) + nn(1:vlen) * H(1:vlen)
    if (id == 0) return

    select case(C_KIND)
    case (C_PBE)
       END = DEC
    case (C_PW91)
       END = DRS
    case default
       STOP
    end select

    do i=DPHI,END

       dHdn(1:vlen,DA) = dHdn(1:vlen,DA) &
            + dH(1:vlen,i) * dXdx(1:vlen,i,DA)
       dHdn(1:vlen,DB) = dHdn(1:vlen,DB) &
            + dH(1:vlen,i) * dXdx(1:vlen,i,DB)
       dHdn(1:vlen,DG) = dHdn(1:vlen,DG) &
            + dH(1:vlen,i) * dXdx(1:vlen,i,DG)

       if (id==2) then
          do j=DPHI,END

             dHdndn(1:vlen,DADA) = dHdndn(1:vlen,DADA) &
                  + d2H(1:vlen,i,j) * dXdx(1:vlen,i,DA) * dXdx(1:vlen,j,DA)
             dHdndn(1:vlen,DADB) = dHdndn(1:vlen,DADB) &
                  + d2H(1:vlen,i,j) * dXdx(1:vlen,i,DA) * dXdx(1:vlen,j,DB)
             dHdndn(1:vlen,DBDB) = dHdndn(1:vlen,DBDB) &
                  + d2H(1:vlen,i,j) * dXdx(1:vlen,i,DB) * dXdx(1:vlen,j,DB)

             dHdndg(1:vlen,DGDA) = dHdndg(1:vlen,DGDA) &
                  + d2H(1:vlen,i,j) * dXdx(1:vlen,i,DG) * dXdx(1:vlen,j,DA)
             dHdndg(1:vlen,DGDB) = dHdndg(1:vlen,DGDB) &
                  + d2H(1:vlen,i,j) * dXdx(1:vlen,i,DG) * dXdx(1:vlen,j,DB)

             dHdgdg(1:vlen) = dHdgdg(1:vlen) &
                  + d2H(1:vlen,i,j) * dXdx(1:vlen,i,DG) * dXdx(1:vlen,j,DG)
          enddo

          dHdndn(1:vlen,DADA) = dHdndn(1:vlen,DADA) &
               + dH(1:vlen,i) * d2Xdx2(1:vlen,i,DA,DA)
          dHdndn(1:vlen,DADB) = dHdndn(1:vlen,DADB) &
               + dH(1:vlen,i) * d2Xdx2(1:vlen,i,DA,DB)
          dHdndn(1:vlen,DBDB) = dHdndn(1:vlen,DBDB) &
               + dH(1:vlen,i) * d2Xdx2(1:vlen,i,DB,DB)

          dHdndg(1:vlen,DGDA) = dHdndg(1:vlen,DGDA) &
               + dH(1:vlen,i) * d2Xdx2(1:vlen,i,DG,DA)
          dHdndg(1:vlen,DGDB) = dHdndg(1:vlen,DGDB) &
               + dH(1:vlen,i) * d2Xdx2(1:vlen,i,DG,DB)

          dHdgdg(1:vlen) = dHdgdg(1:vlen) &
               + dH(1:vlen,i) * d2Xdx2(1:vlen,i,DG,DG)
       end if
    enddo

    if (id==2) then
       dfc_drhodrho(1:vlen,DADA) = dfc_drhodrho(1:vlen,DADA) &
            &                    + nn * dHdndn(1:vlen,DADA) + dHdn(1:vlen,DA) + dHdn(1:vlen,DA)
       dfc_drhodrho(1:vlen,DADB) = dfc_drhodrho(1:vlen,DADB) &
            &                    + nn * dHdndn(1:vlen,DADB) + dHdn(1:vlen,DA) + dHdn(1:vlen,DB)
       dfc_drhodrho(1:vlen,DBDB) = dfc_drhodrho(1:vlen,DBDB) + &
            nn * dHdndn(1:vlen,DBDB) + dHdn(1:vlen,DB) + dHdn(1:vlen,DB)

       dfc_drhodgamma(1:vlen,DGDA) = dfc_drhodgamma(1:vlen,DGDA) &
            &                      + nn * dHdndg(1:vlen,DGDA) + dHdn(1:vlen,DG)
       dfc_drhodgamma(1:vlen,DGDB) = dfc_drhodgamma(1:vlen,DGDB) &
            &                      + nn * dHdndg(1:vlen,DGDB) + dHdn(1:vlen,DG)

       dfc_dgammadgamma(1:vlen)    = dfc_dgammadgamma(1:vlen) &
            &                      + nn * dHdgdg(1:vlen)
    end if

    dfcdn(1:vlen,DA) = dfcdn(1:vlen,DA) + H(1:vlen) + nn(1:vlen) * dHdn(1:vlen,DA)
    dfcdn(1:vlen,DB) = dfcdn(1:vlen,DB) + H(1:vlen) + nn(1:vlen) * dHdn(1:vlen,DB)

    dfcdg(1:vlen   ) = dfcdg(1:vlen   )             + nn(1:vlen) * dHdn(1:vlen,DG)

  end subroutine pbe_or_pw91
  !********************************************************************

  subroutine H_pbe(nd, p, t, e, nv, H, dH, d2H)
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: p,t,e
    real(kind=r8_kind),intent(out),dimension(nv)     :: H
    real(kind=r8_kind),intent(out),dimension(nv,3)   :: dH
    real(kind=r8_kind),intent(out),dimension(nv,3,3) :: d2H
    real   (kind=r8_kind),parameter ::&
         beta  = 0.06672455060314922_r8_kind, &!!0.066725_r8_kind,    &
         gamma = 0.03109069086965489503494086371_r8_kind, &
         bdg   = beta/gamma

    !! DO NOT FORGET t means t^2
    call H0_drv(nd, p, t, e, nv, H, dH, d2H, beta, gamma, bdg)

  end subroutine H_pbe
  !********************************************************************

  subroutine H_pw91(nd, rs, p, t, e, nv, H, dH, d2H)
    use constants, only: one, two, three, four, five, pi
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: rs, p,t,e
    real(kind=r8_kind),intent(out),dimension(nv)     :: H
    real(kind=r8_kind),intent(out),dimension(nv,4)   :: dH
    real(kind=r8_kind),intent(out),dimension(nv,4,4) :: d2H
    real   (kind=r8_kind) :: beta, gamma, bdg

    beta  = nu*Cc0
    gamma = beta * beta / (two * alpha)
    bdg   = beta/gamma

    !! DO NOT FORGET t means t^2
    call H0_drv(nd,     p, t, e, nv, H, dH, d2H, beta, gamma, bdg)
    call H1_drv(nd, rs, p, t,    nv, H, dH, d2H)

  end subroutine H_pw91
  !********************************************************************

  subroutine H0_drv(nd, p, t, e, nv, H, dH, d2H, beta, gamma, bdg)
    use constants
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: p,t,e
    real(kind=r8_kind),intent(out),dimension(nv)     :: H
    real(kind=r8_kind),intent(out),dimension(nv,4)   :: dH
    real(kind=r8_kind),intent(out),dimension(nv,4,4) :: d2H
    real(kind=r8_kind),intent(in) :: beta, gamma, bdg

    real(kind=r8_kind),dimension(nv)             :: p2,p3,Hln,A
    real(kind=r8_kind),dimension(nv,2)           :: dHln, dA
    real(kind=r8_kind),dimension(nv,3)           :: d2Hln,d2A

    p2  = p**two
    p3  = p *p2

    call A_drv  (nd,p,e,nv, A,  dA,  d2A, beta, gamma, bdg  )
    call Hln_drv(nd,t  ,A ,nv, Hln,dHln,d2Hln, beta, gamma, bdg)

    H = gamma*p3*Hln
    if (nd == 0 ) return
    !! 1st derv
    dH(1:nv,DT2  ) = gamma*p3*dHln(1:nv,DT)
    dH(1:nv,DEC  ) = gamma*p3*dHln(1:nv,DX)*dA(1:nv,DE)
    dH(1:nv,DPHI) = gamma*p2*(three*Hln(1:nv)&
         &        + p*dHln(1:nv,DX)*dA(1:nv,DP))
    if (nd == 1) return
    !! 2nd derv
    d2H(1:nv,DT2 ,DT2 ) = gamma*p3*d2Hln(1:nv,DTDT)
    d2H(1:nv,DEC ,DEC ) = gamma*p3*(d2Hln(1:nv,DXDX)*dA (1:nv,DE)**2 &
         &          +           dHln (1:nv,DX  )*d2A(1:nv,DE))
    d2H(1:nv,DPHI,DPHI) = six*gamma*p*  Hln (1:nv) &
         &          + six*gamma*p2 * dHln (1:nv,DX  )*dA (1:nv,DP)    &
         &          + gamma*p3     *(d2Hln(1:nv,DXDX)*dA (1:nv,DP)**2 &
         &          +                dHln (1:nv,DX  )*d2A(1:nv,DPDP))

    d2H(1:nv,DPHI,DT2 ) = gamma*p2*(three*dHln(1:nv,DT) + p*d2Hln(1:nv,DTDX)*dA(1:nv,DP))
    d2H(1:nv,DT2 ,DPHI) = d2H(1:nv,DPHI,DT2)
    d2H(1:nv,DT2 ,DEC ) = gamma*p3*d2Hln(1:nv,DTDX)*dA(1:nv,DE)
    d2H(1:nv,DEC ,DT2 ) = d2H(1:nv,DT2 ,DEC)
    d2H(1:nv,DPHI,DEC ) = gamma*p2*(three*dHln(1:nv,DX)*dA(1:nv,DE)   &
         &          + p*d2Hln(1:nv,DXDX)*dA (1:nv,DE)*dA(1:nv,DP)&
         &          + p*dHln (1:nv,DX  )*d2A(1:nv,DEDP))
    d2H(1:nv,DEC ,DPHI) = d2H(1:nv,DPHI,DEC)

  end subroutine H0_drv
  !********************************************************************

  subroutine H1_drv(nd, rs, p, t, nv, H, dH, d2H)
    use constants
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: rs,p,t
    real(kind=r8_kind),intent(inout),dimension(nv)          :: H
    real(kind=r8_kind),intent(inout),dimension(nv,DPHI:DRS) :: dH
    real(kind=r8_kind),intent(inout),dimension(nv,DPHI:DRS,DPHI:DRS) :: d2H

    real(kind=r8_kind),dimension(nv)                   :: F, S
    real(kind=r8_kind),dimension(nv,DPHI:DRS)          :: dF, dS ![DPHI,DT2,DEC,DRS]
    real(kind=r8_kind),dimension(nv,DPHI:DRS,DPHI:DRS) :: d2F,d2S

    integer(kind=i4_kind) :: i,j


    call F_drv(nd,rs,p,t,nv, F, dF ,d2F)
    call S_drv(nd,rs,p,t,nv, S, dS, d2S)

    H = H + F * S
    if (nd == 0 ) return
    !! 1st derv
    do i = DPHI,DRS
       dH(1:nv,i) = dH(1:nv,i) + dF(1:nv,i)*S(1:nv) + F(1:nv)*dS(1:nv,i)
    end do
    if (nd == 1 ) return
    do i = DPHI,DRS
       do j = DPHI,DRS
          d2H(1:nv,i,j) = d2H(1:nv,i,j) &
               &        + d2F(1:nv,i,j)  *  S(1:nv) &
               &        + dF (1:nv,i)    * dS(1:nv,j) &
               &        + dF (1:nv,j)    * dS(1:nv,i) &
               &        + d2S(1:nv,i,j)  *  F(1:nv)
       end do
    end do

  end subroutine H1_drv
  !********************************************************************

   subroutine F_drv(nd, rs, p, t, nv, F, dF, d2F)
    use constants
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: rs,p,t
    real(kind=r8_kind),intent(out),dimension(nv)          :: F
    real(kind=r8_kind),intent(out),dimension(nv,DPHI:DRS) :: dF
    real(kind=r8_kind),intent(out),dimension(nv,DPHI:DRS,DPHI:DRS) :: d2F

    real(kind=r8_kind),dimension(nv)    :: p3
    real(kind=r8_kind),dimension(nv,3)  :: Cc ![VAL,DRS,DR2]

    call Cc_drv(nd,nv,rs,Cc)

    p3 = p * p * p

    F = nu*(Cc(1:nv,1)-Cc0-three/seven*Cx)*p3*t
    if (nd == 0 ) return
    dF(1:nv,DPHI) = nu*(Cc(1:nv,1)-Cc0-three/seven*Cx) * (three * p * p) * t
    dF(1:nv,DT2 ) = nu*(Cc(1:nv,1)-Cc0-three/seven*Cx) * p3
    dF(1:nv,DRS ) = nu* Cc(1:nv,2) * p3 * t
    dF(1:nv,DEC ) = ZERO
    if (nd == 1 ) return
    d2F(1:nv,DPHI,DPHI) = nu*(Cc(1:nv,1)-Cc0-three/seven*Cx) * (six * p) * t
    d2F(1:nv,DT2 ,DPHI) = nu*(Cc(1:nv,1)-Cc0-three/seven*Cx) * (three * p * p)
    d2F(1:nv,DRS ,DPHI) = nu* Cc(1:nv,2) * t * (three * p * p)
    d2F(1:nv,DEC ,DPHI) = ZERO

    d2F(1:nv,DPHI,DT2 ) = d2F(1:nv,DT2 ,DPHI)
    d2F(1:nv,DT2 ,DT2 ) = ZERO
    d2F(1:nv,DRS ,DT2 ) = nu* Cc(1:nv,2) * p3
    d2F(1:nv,DEC ,DT2 ) = ZERO

    d2F(1:nv,DPHI,DRS ) = d2F(1:nv,DRS ,DPHI)
    d2F(1:nv,DT2 ,DRS ) = d2F(1:nv,DRS ,DT2 )
    d2F(1:nv,DRS ,DRS ) = nu* Cc(1:nv,3) * p3 * t
    d2F(1:nv,DEC ,DRS ) = ZERO

    d2F(1:nv,DPHI,DEC ) = ZERO
    d2F(1:nv,DT2 ,DEC ) = ZERO
    d2F(1:nv,DRS ,DEC ) = ZERO
    d2F(1:nv,DEC ,DEC ) = ZERO

  end subroutine F_drv
  !********************************************************************

  subroutine S_drv(nd, rs, p, t, nv, S, dS, d2S)
    use constants
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: rs,p,t
    real(kind=r8_kind),intent(out),dimension(nv)          :: S
    real(kind=r8_kind),intent(out),dimension(nv,DPHI:DRS) :: dS
    real(kind=r8_kind),intent(out),dimension(nv,DPHI:DRS,DPHI:DRS) :: d2S

    real(kind=r8_kind),dimension(nv)    :: p2, p3,p4, ine
    real(kind=r8_kind)            :: ksf
    real(kind=r8_kind), parameter :: HNDRT = 100.0_r8_kind

    p2 = p  * p
    p3 = p2 * p
    p4 = p3 * p

    ksf = four/(1.91916_r8_kind*pi)
    ine = -HNDRT*p4*t*ksf*rs
    S   = exp(ine)
    if (nd == 0 ) return
    dS(1:nv,DPHI) = four * (-HNDRT*p3*t*ksf*rs) * S
    dS(1:nv,DT2 ) =         -HNDRT*p4*  ksf*rs  * S
    dS(1:nv,DRS ) =         -HNDRT*p4*t*ksf     * S
    dS(1:nv,DEC ) = ZERO
    if (nd == 1 ) return
    d2S(1:nv,DPHI,DPHI) = S * ( three * four * (-HNDRT*p2*t*ksf*rs)  &
         + (four * (-HNDRT*p3*t*ksf*rs)) * (four *(-HNDRT*p3*t*ksf*rs)))
    d2S(1:nv,DT2 ,DPHI) = S * ( four * (-HNDRT*p3*ksf*rs)  &
         + (four * (-HNDRT*p3*t*ksf*rs)) *        (-HNDRT*p4*  ksf*rs) )
    d2S(1:nv,DRS ,DPHI) = S * ( four * (-HNDRT*p3*t*ksf) &
         + (four * (-HNDRT*p3*t*ksf*rs)) *        (-HNDRT*p4*t*ksf   ) )
    d2S(1:nv,DEC ,DPHI) = ZERO

    d2S(1:nv,DPHI,DT2 ) = d2S(1:nv,DT2 ,DPHI)
    d2S(1:nv,DT2 ,DT2 ) = S * (                 (-HNDRT*p4*ksf*rs)**TWO )
    d2S(1:nv,DRS ,DT2 ) = S * (-HNDRT*p4*ksf +  (-HNDRT*p4*ksf*rs) &
         * (-HNDRT*p4*t*ksf) )
    d2S(1:nv,DEC ,DT2 ) = ZERO

    d2S(1:nv,DPHI,DRS ) = d2S(1:nv,DRS ,DPHI)
    d2S(1:nv,DT2 ,DRS ) = d2S(1:nv,DRS ,DT2 )
    d2S(1:nv,DRS ,DRS ) = S * ((-HNDRT*p4*t*ksf)**TWO)
    d2S(1:nv,DEC ,DRS ) = ZERO

    d2S(1:nv,DPHI,DEC ) = ZERO
    d2S(1:nv,DT2 ,DEC ) = ZERO
    d2S(1:nv,DRS ,DEC ) = ZERO
    d2S(1:nv,DEC ,DEC ) = ZERO

  end subroutine S_drv
  !********************************************************************

  subroutine Cc_drv(nd,nv,rs,Cc)
    use constants
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: rs
    real(kind=r8_kind),intent(out)               :: Cc(:,:)

    real(kind=r8_kind),dimension(nv)    :: rs2
    real(kind=r8_kind),dimension(nv)    :: Fup,  Fdn
    real(kind=r8_kind),dimension(nv)    :: dFup, dFdn
    real(kind=r8_kind),dimension(nv)    :: d2Fup, d2Fdn

    rs2 = rs * rs
    Fup = 2.568_r8_kind + 23.266_r8_kind * rs + 0.007389_r8_kind * rs2
    Fdn = one + 8.723_r8_kind * rs + 0.472_r8_kind * rs2 + 0.073890_r8_kind * rs2 *rs
    Cc(:nv,1) = Fup / Fdn / 1000.0_r8_kind - Cx
    if (nd == 0) return
    dFup = 23.266_r8_kind + two * 0.007389_r8_kind * rs
    dFdn =  8.723_r8_kind + two * 0.472_r8_kind    * rs + three * 0.073890_r8_kind * rs2
    Cc(:nv,2) = (dFup*Fdn - dFdn * Fup) / (Fdn*Fdn) / 1000.0_r8_kind
    if (nd == 1) return
    d2Fup = two * 0.007389_r8_kind
    d2Fdn = two * 0.472_r8_kind     + six * 0.073890_r8_kind * rs
    Cc(:nv,3) = (d2Fup * Fdn**two &
         &    -  two * dFup * dFdn * Fdn &
         &    +  two *  Fup * dFdn ** two &
         &    -  Fup * d2Fdn * Fdn ) &
         &    / Fdn ** three

  end subroutine Cc_drv
  !********************************************************************

  subroutine Hln_drv(nd, t, A, nv, H, dH, d2H, beta, gamma, bdg)
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: t,A
    real(kind=r8_kind),intent(out),dimension(nv)   :: H
    real(kind=r8_kind),intent(out),dimension(nv,2) :: dH
    real(kind=r8_kind),intent(out),dimension(nv,3) :: d2H
    real(kind=r8_kind),intent(in)  :: beta, gamma, bdg

    real(kind=r8_kind),dimension(nv)             :: ct, At, At2
    real(kind=r8_kind),dimension(nv)             :: g,dgda,dgdt,dgdtdt,dgdada,dgdadt

    ct  = bdg*t
    At  = one+ A*t
    At2 = At +(A*t)**two
    g      =  ct*(At/At2)
    H = log(one+g)
    if (nd == 0) return
    !! 1st derv
    dgdt   =  bdg*(At+A*t)/At2**two
    dgda   =     -ct*t*t*A*(1+At)/At2**two
    dH(1:nv,DX) = dgda/(1+g)
    dH(1:nv,DT) = dgdt/(1+g)
    if (nd == 1) return
    !! 2nd derv
    dgdadt = -six*ct*t*A*At      /At2**three
    dgdtdt = -six*ct*A*A*At      /At2**three
    dgdada =  two*ct*t*t*(-one+three*(A*t)**two*At)/At2**three
    d2H(1:nv,DXDX) = dgdada/(1+g) - (dgda/(1+g))**two
    d2H(1:nv,DTDT) = dgdtdt/(1+g) - (dgdt/(1+g))**two
    d2H(1:nv,DTDX) = dgdadt/(1+g) - (dgdt*dgda/(1+g)**two)

  end subroutine Hln_drv
  !********************************************************************

  subroutine A_drv(nd, p, e, nv, A, dA, d2A, beta, gamma, bdg)
    integer(kind=i4_kind),intent(in)             :: nd
    integer(kind=i4_kind),intent(in)             :: nv
    real(kind=r8_kind),intent(in),dimension(nv)  :: p,e
    real(kind=r8_kind),intent(out),dimension(nv)   :: A
    real(kind=r8_kind),intent(out),dimension(nv,2) :: dA
    real(kind=r8_kind),intent(out),dimension(nv,3) :: d2A
    real(kind=r8_kind),intent(in)  :: beta, gamma, bdg

    real(kind=r8_kind),dimension(nv)             :: p3, p4, ee, Xe

    p3  = p**three
    p4  = p3*p
    Xe = -e/(gamma*p3)
    ee = exp(Xe)
    A = bdg/((ee-one)+eps)
    if (nd == 0) return
    !! 1st derv
    where (Xe .le. 300.0_r8_kind)
      dA(1:nv,DP) = -three*bdg*e*ee/((((ee-one)**two)*p4*gamma)+eps)
      dA(1:nv,DE) =        bdg*ee  /((((ee-one)**two)*p3*gamma)+eps)
    elsewhere
      dA(1:nv,DP) = 0.0_rk
      dA(1:nv,DE) = 0.0_rk
    endwhere
    if (nd == 1) return
    !! 2nd derv
    where (Xe .le. 300.0_r8_kind)
      d2A(1:nv,DPDP) = three*bdg*e*ee*(three*e*ee+four*gamma*p3*ee-four*gamma*p3+three*e)&
            &           /(((ee-1)**three*(p4*gamma)**two)+eps)
      d2A(1:nv,DEDE) = bdg*ee*(ee+1)/(((ee-1)**three*(gamma*p3)**two)+eps)
      d2A(1:nv,DEDP) = - three*bdg*ee*(e*ee+gamma*p3*ee-gamma*p3+e)&
            &           /(((ee-1)**three*gamma**two*p**seven)+eps)
    elsewhere
      d2A(1:nv,DPDP) = 0.0_rk
      d2A(1:nv,DEDE) = 0.0_rk
      d2A(1:nv,DEDP) = 0.0_rk
    endwhere


  end subroutine A_drv
  !********************************************************************

  !--------------- End of module -------------------------------------
end module gga_response_module
